//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include "image.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <petscdmda.h>
#include <petscvec.h>

#include "fd_routines.hpp"
#include "indexing.hpp"
#include "iterator_routines.hpp"
#include "map.hpp"
#include "math_utils.hpp"
#include "exceptions.hpp"

#include "baseloader.hpp"

integer ImageBase::instance_id_counter = 0;

// Public Methods

ImageBase::ImageBase(const intvector& shape, MPI_Comm comm)
  : m_comm(comm), m_ndim(shape.size()),
    m_shape(shape), // const on shape causes copy assignment (c++11)
    m_localvec(create_unique_vec()), m_globalvec(create_unique_vec()), m_dmda(create_shared_dm()),
    instance_id(instance_id_counter++)
{
  if (m_shape.size() != 3)
  {
    if (m_shape.size() == 2)
    {
      m_shape.push_back(1);
    }
    else
    {
      throw InternalError("image shape should be 2D or 3D", __FILE__, __LINE__);
    }
  }
  if (m_shape[2] == 1)
  {
    m_ndim = 2;
  }

  initialize_dmda();
  initialize_vectors();
}

const floating* ImageBase::get_raw_data_ro() const
{
  const floating* ptr;
  PetscErrorCode perr = VecGetArrayRead(*m_globalvec, &ptr);
  CHKERRXX(perr);
  return ptr;
}

void ImageBase::release_raw_data_ro(const floating*& ptr) const
{
  PetscErrorCode perr = VecRestoreArrayRead(*m_globalvec, &ptr);
  CHKERRXX(perr);
}

Vec_unique ImageBase::get_raw_data_row_major() const
{
  Vec_unique rmlocalpart = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*m_globalvec, rmlocalpart.get());
  CHKERRXX(perr);

  intvector locs(3, 0);
  intvector widths(3, 0);
  perr = DMDAGetCorners(*m_dmda, &locs[0], &locs[1], &locs[2], &widths[0], &widths[1], &widths[2]);
  CHKERRXX(perr);

  integer ownedlo;
  perr = VecGetOwnershipRange(*m_globalvec, &ownedlo, nullptr);
  integer localsize = std::accumulate(widths.begin(), widths.end(), 1, std::multiplies<>());
  intvector cmidxn(localsize);
  std::iota(cmidxn.begin(), cmidxn.end(), 0);

  std::transform(
      cmidxn.begin(), cmidxn.end(), cmidxn.begin(), [widths, ownedlo](integer idx) -> integer {
        return idx_cmaj_to_rmaj(idx, widths) + ownedlo;
      });

  IS_unique src_is = create_unique_is();
  perr = ISCreateStride(m_comm, localsize, ownedlo, 1, src_is.get());
  CHKERRXX(perr);
  IS_unique tgt_is = create_unique_is();
  perr = ISCreateGeneral(m_comm, localsize, cmidxn.data(), PETSC_USE_POINTER, tgt_is.get());
  CHKERRXX(perr);
  VecScatter_unique sct = create_unique_vecscatter();
  perr = VecScatterCreate(*m_globalvec, *src_is, *rmlocalpart, *tgt_is, sct.get());
  CHKERRXX(perr);

  perr = VecScatterBegin(*sct, *m_globalvec, *rmlocalpart, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);
  perr = VecScatterEnd(*sct, *m_globalvec, *rmlocalpart, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);

  return rmlocalpart;
}

floating ImageBase::masked_normalize(const Mask& mask)
{
  Vec_unique tmp = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*m_globalvec, tmp.get());CHKERRXX(perr);
  perr = VecPointwiseMult(*tmp, *mask.global_vec(), *m_globalvec);CHKERRXX(perr);

  floating norm;
  perr = VecSum(*tmp, &norm);CHKERRXX(perr);
  norm = mask.npoints() / norm;
  perr = VecScale(*m_globalvec, norm);CHKERRXX(perr);

  return norm;
}

void ImageBase::copy_data(const ImageBase &img)
{
  PetscErrorCode perr = VecCopy(*img.global_vec(), *m_globalvec);
  CHKERRXX(perr);
  update_local_from_global();
}

Vec_unique ImageBase::gradient(integer dim)
{
  // New global vec must be a duplicate of image global
  Vec_shared grad = create_shared_vec();
  PetscErrorCode perr = VecDuplicate(*(m_globalvec), grad.get());
  CHKERRXX(perr);

  // Ensure we have up to date ghost cells
  DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  // Nothing really useful to interleave
  DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);

  return fd::gradient_to_global_unique(*m_dmda, *m_localvec, dim);
}

Vec_unique ImageBase::scatter_to_zero() const
{
  Vec_unique new_vec = create_unique_vec();
  VecScatter_unique sct = create_unique_vecscatter();
  VecScatterCreateToZero(*m_globalvec, sct.get(), new_vec.get());
  VecScatterBegin(*sct, *m_globalvec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(*sct, *m_globalvec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);

  return new_vec;
}

// Return scale factor
floating ImageBase::normalize()
{
  floating norm;
  PetscErrorCode perr = VecMax(*m_globalvec, nullptr, &norm);
  CHKERRXX(perr);
  norm = 1.0 / norm;
  perr = VecScale(*m_globalvec, norm);
  CHKERRXX(perr);
  return norm;
}

void ImageBase::binarize()
{
  binarize_vector(*m_globalvec, 0.01);
  update_local_from_global();
}

//// Protected Methods

ImageBase::ImageBase(const ImageBase& imb)
  : m_comm(imb.m_comm), m_ndim(imb.m_ndim), m_shape(imb.m_shape),
    m_localvec(create_shared_vec()), m_globalvec(create_shared_vec()), m_dmda(imb.m_dmda),
    instance_id(instance_id_counter++)
{
  initialize_vectors();
}

void ImageBase::initialize_dmda()
{
  integer dof_per_node = 1;
  integer stencil_width = 1;
  PetscErrorCode perr;

  // Make sure things get gracefully cleaned up
  m_dmda = create_shared_dm();
  perr = DMDACreate3d(m_comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, // BCs
      DMDA_STENCIL_STAR,                        // stencil shape
      m_shape[0], m_shape[1], m_shape[2],       // global mesh shape
      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, // ranks per dim
      dof_per_node, stencil_width,              // dof per node, stencil size
      nullptr, nullptr, nullptr,                // partition sizes nullptr -> petsc chooses
      m_dmda.get());
  CHKERRXX(perr);

  perr = DMSetUp(*(m_dmda));
  CHKERRXX(perr);
}

void ImageBase::initialize_vectors()
{
  PetscErrorCode perr;
  m_localvec = create_shared_vec();
  m_globalvec = create_shared_vec();
  perr = DMCreateGlobalVector(*m_dmda, m_globalvec.get());
  CHKERRXX(perr);
  debug_creation(*m_globalvec, std::string("image_global_") + std::to_string(instance_id));
  perr = DMCreateLocalVector(*m_dmda, m_localvec.get());
  CHKERRXX(perr);
  debug_creation(*m_localvec, std::string("image_local_") + std::to_string(instance_id));
}

void ImageBase::update_local_from_global()
{
  PetscErrorCode perr = DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  CHKERRXX(perr);
  perr = DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  CHKERRXX(perr);
}

Vec_unique ImageBase::get_raw_data_natural() const
{
  Vec_unique natvector = create_unique_vec();
  PetscErrorCode perr = DMDACreateNaturalVector(*m_dmda, natvector.get());
  CHKERRXX(perr);

  perr = DMDAGlobalToNaturalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *natvector);
  CHKERRXX(perr);
  perr = DMDAGlobalToNaturalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *natvector);
  CHKERRXX(perr);

  return natvector;
}

floating ImageBase::mutual_information(const ImageBase &other)
{
  // Mutual information is given by integrating P(X, Y)log(P(X,Y)/(P(X)P(Y)) over all X and all Y,
  // In this case the probabilities can be calculated from a 2D histogram of X against Y.

  // Could also just find 2D and perform summations after comm, but for now do it this way
  floatvector xhist(mi_resolution, 0);
  floatvector yhist(mi_resolution, 0);
  floatvector2d xyhist(mi_resolution, floatvector(mi_resolution, 0));

  integer img_localsize;
  PetscErrorCode perr = VecGetLocalSize(*m_globalvec, &img_localsize);
  CHKERRXX(perr);

  floating max1, max2;
  perr = VecMax(*m_globalvec, nullptr, &max1);CHKERRXX(perr);
  perr = VecMax(*other.global_vec(), nullptr, &max2);CHKERRXX(perr);
  floating max = std::max(max1, max2);

  // Only need RO data from PETSc vecs
  floating const *x_data, *y_data;
  perr = VecGetArrayRead(*m_globalvec, &x_data);
  CHKERRXX(perr);
  perr = VecGetArrayRead(*other.global_vec(), &y_data);
  CHKERRXX(perr);

  for (integer idx = 0; idx < img_localsize; idx++)
  {
    integer x = std::lround(x_data[idx]/max * static_cast<double>(mi_resolution-1));
    integer y = std::lround(y_data[idx]/max * static_cast<double>(mi_resolution-1));
    xhist[x] += 1;
    yhist[y] += 1;
    xyhist[x][y] += 1;
  }

  int rank;
  MPI_Comm_rank(m_comm, &rank);
  if(rank == 0)
  {
    MPI_Reduce(MPI_IN_PLACE, xhist.data(), xhist.size(), MPIU_SCALAR, MPI_SUM, 0, m_comm);
    MPI_Reduce(MPI_IN_PLACE, yhist.data(), yhist.size(), MPIU_SCALAR, MPI_SUM, 0, m_comm);
    MPI_Reduce(MPI_IN_PLACE, xyhist.data(), xyhist.size(), MPIU_SCALAR, MPI_SUM, 0, m_comm);
  } 
  else
  {
    MPI_Reduce(xhist.data(), xhist.data(), xhist.size(), MPIU_SCALAR, MPI_SUM, 0, m_comm);
    MPI_Reduce(yhist.data(), yhist.data(), yhist.size(), MPIU_SCALAR, MPI_SUM, 0, m_comm);
    MPI_Reduce(xyhist.data(), xyhist.data(), xyhist.size(), MPIU_SCALAR, MPI_SUM, 0, m_comm);
  }

  // Need all probability distributions to sum to 1.0
  integer pix_tot = this->size();
  std::transform(xhist.cbegin(), xhist.cend(), xhist.begin(),
                 [pix_tot](floating x) -> floating {return x/pix_tot;});
  std::transform(xhist.cbegin(), xhist.cend(), xhist.begin(),
                 [pix_tot](floating y) -> floating {return y/pix_tot;});
  for(auto &vec1d: xyhist)
  {
    std::transform(vec1d.cbegin(), vec1d.cend(), vec1d.begin(),
                   [pix_tot](floating y) -> floating {return y/pix_tot;});
  }

  // Use Kahan summation to try to minimize error
  floating mi_total(0);
  floating mi_err(0);
  if(rank == 0)
  {
    for(size_t ix(0); ix < xhist.size(); ix++)
    {
      for(size_t iy(0); iy < yhist.size(); iy++)
      {
        if(xyhist[ix][iy] > 0)
        {
          floating mi_part = xyhist[ix][iy] * std::log(xyhist[ix][iy]/(xhist[ix]*yhist[iy]));
          mi_part -= mi_err;
          floating tmp = mi_total + mi_part;
          mi_err = (tmp - mi_total) - mi_part;
          mi_total = tmp;
        }
      }
    }
  }

  MPI_Bcast(&mi_total, 1, MPIU_SCALAR, 0, m_comm);

  return mi_total;
}
