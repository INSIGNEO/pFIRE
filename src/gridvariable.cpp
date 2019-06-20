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

#include "exceptions.hpp"
#include "gridvariable.hpp"

#include "indexing.hpp"

GridVariable::GridVariable(const GridVariable &pattern)
  : _comm(pattern.comm()), _commsize(pattern.commsize()), _shape(pattern.shape()),
    _index_max(_shape - 1), _ndim(pattern.ndim()), _ndof(pattern.ndof()), _dmda(pattern.dmda_ptr())
{
  initialise_vectors();
  initialise_additional_info();
}

GridVariable::GridVariable(const intcoord &shape, const integer &ndof, const MPI_Comm &comm)
  : _comm(comm), _shape(shape), _index_max(_shape - 1),
    _ndim(shape.back() < 2 ? shape.size() - 1 : shape.size()), _ndof(ndof)
{
  MPI_Comm_size(_comm, &_commsize);

  initialise_dmda();
  initialise_vectors();
  initialise_additional_info();
}

std::pair<intcoord, intcoord> GridVariable::owned_range() const
{
  return std::make_pair(_local_offset, _local_offset + _local_shape);
}

Vec_unique GridVariable::scatter_to_zero() const
{
  Vec_unique new_vec = create_unique_vec();
  VecScatter_unique sct = create_unique_vecscatter();
  VecScatterCreateToZero(*_globalvec, sct.get(), new_vec.get());
  VecScatterBegin(*sct, *_globalvec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(*sct, *_globalvec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);

  return new_vec;
}

void GridVariable::initialise_dmda()
{
  // clang-format off
  _dmda = create_shared_dm();
  PetscErrorCode perr = DMDACreate3d(_comm,
                                     _bnd_type, _bnd_type, _bnd_type,          // these are globally sensible
                                     _stencil_type,                            // and here too
                                     _shape[0], _shape[1], _shape[2],          // set shape
                                     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, // rank splitting
                                     _ndof, _stencil_width,                    // self explanatory
                                     nullptr, nullptr, nullptr,                // petsc chooses nodes per rank
                                     _dmda.get());
  CHKERRXX(perr);

  perr = DMSetUp(*_dmda);
  CHKERRXX(perr);
}

floatpair GridVariable::minmax() const
{
  floating min, max;
  PetscErrorCode perr = VecMin(this->global_vector(), nullptr, &min);
  CHKERRXX(perr);
  perr = VecMax(this->global_vector(), nullptr, &max);
  CHKERRXX(perr);

  return std::make_pair(min, max);
}

void GridVariable::initialise_additional_info()
{
  // Now make a note of edge locations for later
  PetscErrorCode perr = DMDAGetInfo(*_dmda,
                                    nullptr, nullptr, nullptr, nullptr, // just want the ranks-per-dim arrays
                                    &_ranks_per_dim[0], &_ranks_per_dim[1], &_ranks_per_dim[2],
                                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  CHKERRXX(perr);
  // clang-format on

  _rank_edges = {intvector(_ranks_per_dim[0] + 1, 0), intvector(_ranks_per_dim[1] + 1, 0),
      intvector(_ranks_per_dim[2] + 1, 0)};
  coord<const integer *> chunksizes;
  // The following is not a memory leak, we are picking up a pointer to PETSc internal
  // structures...
  perr = DMDAGetOwnershipRanges(*_dmda, &chunksizes[0], &chunksizes[1], &chunksizes[2]);
  CHKERRXX(perr);

  for (int i = 0; i < 3; i++)
  {
    // Pointers are iterators too (sort of!)
    std::partial_sum(chunksizes[i], chunksizes[i] + _ranks_per_dim[i], _rank_edges[i].begin() + 1);
  }

  DMDALocalInfo info;
  DMDAGetLocalInfo(*_dmda, &info);
  _local_shape = {info.xm, info.ym, info.zm};
  _local_offset = {info.xs, info.ys, info.zs};
}

void GridVariable::initialise_vectors()
{
  _globalvec = create_unique_vec();
  PetscErrorCode perr = DMCreateGlobalVector(*_dmda, _globalvec.get());
  CHKERRXX(perr);
  perr = VecSet(*_globalvec, 0.);
  CHKERRXX(perr);

  _localvec = create_unique_vec();
  perr = DMCreateLocalVector(*_dmda, _localvec.get());
  CHKERRXX(perr);
  perr = VecSet(*_localvec, 0.);
  CHKERRXX(perr);
}

void GridVariable::copy_data_from(const GridVariable &src)
{
  if (this->dmda() != src.dmda())
  {
    throw InternalError("GridVariables must share a DMDA", __FILE__, __LINE__);
  }

  // Do the vector copy
  PetscErrorCode perr = VecCopy(src.global_vector(), this->global_vector());
  CHKERRXX(perr);

  this->update_local_vector();
}

void GridVariable::update_local_vector() const
{
  PetscErrorCode perr = DMGlobalToLocalBegin(
      this->dmda(), this->global_vector(), INSERT_VALUES, this->local_vector());
  CHKERRXX(perr);
  perr =
      DMGlobalToLocalEnd(this->dmda(), this->global_vector(), INSERT_VALUES, this->local_vector());
  CHKERRXX(perr);
}

intcoord GridVariable::find_rank_distribution(const intcoord &meshshape, const integer &commsize)
{
  const integer &M = meshshape[0];
  const integer &N = meshshape[1];
  const integer &P = meshshape[2];

  intcoord rank_distribution;
  integer &m = rank_distribution[0];
  integer &n = rank_distribution[1];
  integer &p = rank_distribution[2];

  integer pm; // p*m
  // This is the same algorithm used in PETSc (as of 3.13)
  // Try to choose squarish blocks on each rank
  n = std::lround(std::pow(commsize * static_cast<floating>(N * N) / (P * M), 1. / 3));
  if (n < 1)
  {
    n = 1;
  }

  while (n > 0)
  {
    pm = commsize / n;
    if (n * pm == commsize)
    {
      break;
    }
    n--;
  }

  m = std::lround(std::sqrt(pm * static_cast<floating>(M) / P));
  if (m < 1)
  {
    m = 1;
  }

  while (m > 0)
  {
    p = pm / m;
    if (p * m == pm)
    {
      break;
    }
    m--;
  }
  if (M > P && m < p)
  {
    integer _m = m;
    m = p;
    p = _m;
  }

  return rank_distribution;
}
