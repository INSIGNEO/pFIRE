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
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <petscdmda.h>
#include <petscvec.h>

#include "basis.hpp"
#include "fd_routines.hpp"
#include "indexing.hpp"
#include "infix_iterator.hpp"
#include "iterator_routines.hpp"
#include "map.hpp"

#include "baseloader.hpp"

integer Image::instance_id_counter = 0;

// Public Methods

Image::Image(const intvector& shape, MPI_Comm comm)
  : m_comm(comm), m_ndim(shape.size()),
    m_shape(shape), // const on shape causes copy assignment (c++11)
    m_localvec(create_unique_vec()), m_globalvec(create_unique_vec()), m_dmda(create_shared_dm()),
    instance_id(instance_id_counter++)
{
  MPI_Comm_rank(comm, &rank);
  if (m_shape.size() != 3)
  {
    if (m_shape.size() == 2)
    {
      m_shape.push_back(1);
    }
    else
    {
      throw std::runtime_error("image shape should be 2D or 3D");
    }
  }
  if (m_shape[2] == 1)
  {
    m_ndim = 2;
  }

  initialize_dmda();
  initialize_vectors();
}

std::unique_ptr<Image> Image::duplicate() const
{
  return std::unique_ptr<Image>(new Image(*this));
}

std::unique_ptr<Image> Image::copy() const
{
  // private copy c'tor prohibits use of std::make_unique, would otherwise do:
  // std::unique_ptr<Image> new_img = std::make_unique<Image>(*this);
  std::unique_ptr<Image> new_img(new Image(*this));

  PetscErrorCode perr = VecCopy(*m_localvec, *new_img->m_localvec);
  CHKERRABORT(m_comm, perr);
  perr = VecCopy(*m_globalvec, *new_img->m_globalvec);
  CHKERRABORT(m_comm, perr);

  return new_img;
}

std::unique_ptr<Image> Image::load_file(
    const std::string& path, const Image* existing, MPI_Comm comm)
{
  BaseLoader_unique loader = BaseLoader::find_loader(path, comm);

  // if image passed assert sizes match and duplicate, otherwise create new image given size
  std::unique_ptr<Image> new_image;
  if (existing != nullptr)
  {
    comm = existing->comm();
    if (!all_true(loader->shape().begin(), loader->shape().end(), existing->shape().begin(),
            existing->shape().end(), std::equal_to<>()))
    {
      throw std::runtime_error("New image must have same shape as existing");
    }
    new_image = existing->duplicate();
  }
  else
  {
    new_image = std::make_unique<Image>(loader->shape(), comm);
  }

  intvector shape(3, 0), offset(3, 0);
  PetscErrorCode perr = DMDAGetCorners(
      *new_image->dmda(), &offset[0], &offset[1], &offset[2], &shape[0], &shape[1], &shape[2]);
  CHKERRABORT(comm, perr);
  // std::transform(shape.cbegin(), shape.cend(), offset.cbegin(), shape.begin(), std::minus<>());

  floating*** vecptr(nullptr);
  perr = DMDAVecGetArray(*new_image->dmda(), *new_image->global_vec(), &vecptr);
  CHKERRABORT(comm, perr);
  loader->copy_scaled_chunk(vecptr, shape, offset);
  perr = DMDAVecRestoreArray(*new_image->dmda(), *new_image->global_vec(), &vecptr);
  CHKERRABORT(comm, perr);

  return new_image;
}
// Return scale factor
floating Image::normalize()
{
  floating norm;
  PetscErrorCode perr = VecSum(*m_globalvec, &norm);
  CHKERRABORT(m_comm, perr);
  norm = this->size() / norm;
  perr = VecScale(*m_globalvec, norm);
  CHKERRABORT(m_comm, perr);
  return norm;
}

// Protected Methods

Image::Image(const Image& image)
  : m_comm(image.m_comm), m_ndim(image.m_ndim), m_shape(image.m_shape),
    m_localvec(create_shared_vec()), m_globalvec(create_shared_vec()), m_dmda(image.m_dmda),
    instance_id(instance_id_counter++)
{
  initialize_vectors();
}

void Image::initialize_dmda()
{
  // Should never change
  DMBoundaryType btype = DM_BOUNDARY_MIRROR; // mirror with 1 point means zero gradient at edges
  DMDAStencilType stype = DMDA_STENCIL_STAR;
  integer stencil_width = 1;
  integer dof_per_node = 1;

  // Calculate rank distribution
  int comm_size, rank;
  MPI_Comm_size(m_comm, &comm_size);
  MPI_Comm_rank(m_comm, &rank);
  auto splitinfo = find_proc_split(m_shape, comm_size);
  ranks_per_dim = splitinfo.first;
  strides = splitinfo.second;

  if (rank == 0)
  {
    std::cout << "Found proc split (";
    std::copy(ranks_per_dim.cbegin(), ranks_per_dim.cend(),
        infix_ostream_iterator<integer>(std::cout, ", "));
    std::cout << ") with strides (";
    std::copy(strides.cbegin(), strides.cend(), infix_ostream_iterator<integer>(std::cout, ", "));
    std::cout << ").\n";
  }

  // Calculate nodes per cell explicitly
  intvector2d nodes_per_cell(0);
  nodes_per_cell.reserve(3);
  // Construct nodes per cell info for DMDACreate3D
  for (integer ridx = 0; ridx < 3; ridx++)
  {
    intvector nodes(ranks_per_dim[ridx], 0);
    std::fill(nodes.begin(), nodes.end(), strides[ridx]);
    *nodes.rbegin() +=
        m_shape[ridx]
        - (ranks_per_dim[ridx] * strides[ridx]); // last cell contains small remainder
    nodes_per_cell.emplace_back(nodes);
  }

  // use a shared ptr to create dmda
  m_dmda = create_shared_dm();
  PetscErrorCode perr;
  perr = DMDACreate3d(m_comm, btype, btype, btype, stype,   // BCs, stencil shape
      m_shape[0], m_shape[1], m_shape[2],                   // global mesh shape
      ranks_per_dim[0], ranks_per_dim[1], ranks_per_dim[2], // ranks per dim
      dof_per_node, stencil_width,                          // dof per node, stencil size
      nodes_per_cell[0].data(), nodes_per_cell[1].data(),
      nodes_per_cell[2].data(), // partition sizes
      m_dmda.get());
  CHKERRABORT(m_comm, perr);

  perr = DMSetUp(*(m_dmda));
  CHKERRABORT(m_comm, perr);

  // finally, map edges to ranks:
  intvector rankloc(3, 0);
  perr = DMDAGetCorners(*m_dmda, &rankloc[0], &rankloc[1], &rankloc[2], nullptr, nullptr, nullptr);
  CHKERRABORT(m_comm, perr);

  rankloc[0] = rankloc[0] / strides[0];
  rankloc[1] = rankloc[1] / strides[1];
  rankloc[2] = rankloc[2] / strides[2];

  // TODO: improve me
  rankmapping.resize(comm_size);
  std::fill(rankmapping.begin(), rankmapping.end(), 0);
  intvector tmpmapping(comm_size, 0);
  PetscSynchronizedPrintf(m_comm, "Rank: %d - Loc: %d\n", rank, ravel(rankloc, ranks_per_dim));
  PetscSynchronizedFlush(m_comm, PETSC_STDOUT);
  tmpmapping[ravel(rankloc, ranks_per_dim)] = rank;

  MPI_Allreduce(
      tmpmapping.data(), rankmapping.data(), rankmapping.size(), MPIU_INT, MPI_SUM, m_comm);
  if (rank == 0)
  {
    std::copy(
        rankmapping.cbegin(), rankmapping.cend(), std::ostream_iterator<integer>(std::cout, " "));
  }
}

/**
 * Return mpi rank of process with region containing "loc"
 * Locations outside of the image bounds return the rank containing the closest valid pixel
 *
 * @param loc location with in image (pixel units)
 *
 * @return rank of process containing "loc"
 */
integer Image::get_rank_of_loc(const floatloc& loc) const
{
  intvector intloc(3, 0);
  n_ary_transform([](floating coord, integer stride) -> integer { return coord / stride; },
      intloc.begin(), loc.cbegin(), loc.cend(), strides.cbegin());

  // Clamp out of bounds values to edge ranks
  for (size_t idx = 0; idx < loc.size(); idx++)
  {
    if (intloc[idx] < 0)
    {
      intloc[idx] = 0;
    }
    else if (intloc[idx] >= ranks_per_dim[idx])
    {
      intloc[idx] = ranks_per_dim[idx] - 1;
    }
  }

  return ravel(intloc, ranks_per_dim);
}

void Image::initialize_vectors()
{
  PetscErrorCode perr;
  m_localvec = create_shared_vec();
  m_globalvec = create_shared_vec();
  perr = DMCreateGlobalVector(*m_dmda, m_globalvec.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*m_globalvec, std::string("image_global_") + std::to_string(instance_id));
  perr = DMCreateLocalVector(*m_dmda, m_localvec.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*m_localvec, std::string("image_local_") + std::to_string(instance_id));
}

void Image::update_local_from_global() const
{
  PetscErrorCode perr = DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
}

Vec_unique Image::gradient(integer dim)
{
  // New global vec must be a duplicate of image global
  Vec_shared grad = create_shared_vec();
  PetscErrorCode perr = VecDuplicate(*(m_globalvec), grad.get());
  CHKERRABORT(m_comm, perr);

  // Ensure we have up to date ghost cells
  DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  // Nothing really useful to interleave
  DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);

  return fd::gradient_to_global_unique(*m_dmda, *m_localvec, dim);
}

Vec_unique Image::scatter_to_zero() const
{
  Vec_unique new_vec = create_unique_vec();
  VecScatter_unique sct = create_unique_vecscatter();
  VecScatterCreateToZero(*m_globalvec, sct.get(), new_vec.get());
  VecScatterBegin(*sct, *m_globalvec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(*sct, *m_globalvec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);

  return new_vec;
}

const floating* Image::get_raw_data_ro() const
{
  const floating* ptr;
  PetscErrorCode perr = VecGetArrayRead(*m_globalvec, &ptr);
  CHKERRABORT(m_comm, perr);
  return ptr;
}

void Image::release_raw_data_ro(const floating*& ptr) const
{
  PetscErrorCode perr = VecRestoreArrayRead(*m_globalvec, &ptr);
  CHKERRABORT(m_comm, perr);
}

Vec_unique Image::get_raw_data_natural() const
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

Vec_unique Image::get_raw_data_row_major() const
{
  Vec_unique rmlocalpart = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*m_globalvec, rmlocalpart.get());
  CHKERRABORT(m_comm, perr);

  intvector locs(3, 0);
  intvector widths(3, 0);
  perr = DMDAGetCorners(*m_dmda, &locs[0], &locs[1], &locs[2], &widths[0], &widths[1], &widths[2]);
  CHKERRABORT(m_comm, perr);

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
  CHKERRABORT(m_comm, perr);
  IS_unique tgt_is = create_unique_is();
  perr = ISCreateGeneral(m_comm, localsize, cmidxn.data(), PETSC_USE_POINTER, tgt_is.get());
  CHKERRABORT(m_comm, perr);
  VecScatter_unique sct = create_unique_vecscatter();
  perr = VecScatterCreate(*m_globalvec, *src_is, *rmlocalpart, *tgt_is, sct.get());
  CHKERRABORT(m_comm, perr);

  perr = VecScatterBegin(*sct, *m_globalvec, *rmlocalpart, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(m_comm, perr);
  perr = VecScatterEnd(*sct, *m_globalvec, *rmlocalpart, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(m_comm, perr);

  return rmlocalpart;
}

ImageInterpolator::ImageInterpolator(const Image& image) : image(image)
{
  image.update_local_from_global();
  PetscErrorCode perr = DMDAVecGetArrayRead(*image.dmda(), *image.local_vec(), &localdata);
  CHKERRXX(perr);
}

ImageInterpolator::~ImageInterpolator()
{
  PetscErrorCode perr = DMDAVecRestoreArrayRead(*image.dmda(), *image.local_vec(), &localdata);
  CHKERRXX(perr);
}

floating ImageInterpolator::operator()(floatloc loc) const
{
  // Check if image is 2D and limit to 2D interpolation if so
  integer npoints = 8;
  if (image.ndim() == 2)
  {
    npoints = 4;
  }

  n_ary_transform([](floating x, integer w) -> floating { return clamp_to_edge(x, w); },
      loc.begin(), loc.cbegin(), loc.cend(), image.shape().cbegin());
  // find low corner
  intloc loc_floor;
  std::transform(loc.cbegin(), loc.cend(), loc_floor.begin(),
      [](floating x) -> integer { return static_cast<integer>(std::floor(x)); });

  // iterate over all corners and add coeff*val to result
  floating result = 0;
  for (integer ipoint = 0; ipoint < npoints; ipoint++)
  {
    intloc ploc = loc_floor;
    for (uinteger idim = 0; idim < image.ndim(); idim++)
    {
      if ((1 << idim) & ipoint)
      {
        ploc[idim] += 1;
      }
    }
    floating coeff = calculate_basis_coefficient(loc.cbegin(), loc.cend(), ploc.cbegin());
    result += coeff * localdata[ploc[2]][ploc[1]][ploc[0]];
  }

  return result;
}
