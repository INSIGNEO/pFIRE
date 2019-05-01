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

#include "petsc_helpers.hpp"

#include <algorithm>
#include <sstream>
#include <fenv.h>

#include "exceptions.hpp"

floating diagonal_sum(const Mat& matrix)
{
  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);

  Vec_unique diag = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(matrix, diag.get(), nullptr);
  CHKERRXX(perr);

  perr = MatGetDiagonal(matrix, *diag);
  CHKERRXX(perr);

  floating diagsum;
  perr = VecSum(*diag, &diagsum);

  return diagsum;
}

bool vecs_equivalent(const Vec &vec1, const Vec &vec2)
{
  // Ensure two vectors have the same local size, global size. and comm.
  MPI_Comm comm1, comm2;
  PetscErrorCode perr = PetscObjectGetComm((PetscObject)vec1, &comm1);
  CHKERRXX(perr);
  perr = PetscObjectGetComm((PetscObject)vec2, &comm2);
  CHKERRXX(perr);
  if (comm1 != comm2)
  {
    return false;
  }

  integer globalsize1, globalsize2;
  perr = VecGetSize(vec1, &globalsize1);
  CHKERRXX(perr);
  perr = VecGetSize(vec2, &globalsize2);
  CHKERRXX(perr);
  if (globalsize1 != globalsize2)
  {
    return false;
  }

  integer localsize1, localsize2;
  perr = VecGetLocalSize(vec1, &localsize1);
  CHKERRXX(perr);
  perr = VecGetLocalSize(vec2, &localsize2);
  CHKERRXX(perr);
  if (localsize1 != localsize2)
  {
    return false;
  }

  return true;
}

void block_precondition(const Mat& normmat, integer size, integer ndim)
{
  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(normmat), &comm);
  // Normalize luminance block of matrix to spatial blocks using diagonal norm
  Vec_unique diag = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(normmat, diag.get(), nullptr);
  CHKERRXX(perr);
  debug_creation(*diag, "Vec_block_normalise");
  perr = MatGetDiagonal(normmat, *diag);
  CHKERRXX(perr);

  // index where spatial dims stop and intensity dim starts
  integer crit_idx = size * ndim;

  // Find rank-local sums for spatial and luminance blocks
  integer rowstart, rowend;
  perr = VecGetOwnershipRange(*diag, &rowstart, &rowend);
  CHKERRXX(perr);
  integer localsize = rowend - rowstart;
  floatvector norm(2, 0.0); // norm[0] is spatial, norm[1] is luminance
  integer spt_end = crit_idx - rowstart;
  spt_end = (spt_end < 0) ? 0 : (spt_end > localsize) ? localsize : spt_end;

  floating* ptr;
  perr = VecGetArray(*diag, &ptr);
  CHKERRXX(perr);
  for (integer idx = 0; idx < spt_end; idx++)
  {
    norm[0] += ptr[idx];
  }
  for (integer idx = spt_end; idx < localsize; idx++)
  {
    norm[1] += ptr[idx];
  }
  perr = VecRestoreArray(*diag, &ptr);
  CHKERRXX(perr);

  // MPI_AllReduce to sum over all processes
  MPI_Allreduce(MPI_IN_PLACE, norm.data(), 2, MPI_DOUBLE, MPI_SUM, comm);

  // calculate average of norms and scaling factor
  norm[0] /= crit_idx;
  norm[1] /= size;
  floating lum_scale = norm[0] / norm[1];

  // reuse diag vector to hold scaling values
  perr = VecGetArray(*diag, &ptr);
  CHKERRXX(perr);
  for (integer idx = 0; idx < spt_end; idx++)
  {
    ptr[idx] = 1;
  }
  for (integer idx = spt_end; idx < localsize; idx++)
  {
    ptr[idx] = lum_scale;
  }
  perr = VecRestoreArray(*diag, &ptr);
  CHKERRXX(perr);

  perr = MatDiagonalScale(normmat, *diag, nullptr);
  CHKERRXX(perr);
}

floating get_condnum_by_poweriter(const Mat& matrix, floating conv_thres, integer max_iter)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO);
  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);
  floating eigen_lo = get_eigenvalue_by_poweriter(matrix, conv_thres, max_iter);

  // Get smallest eigenvalue by spectral shift
  Mat_unique e_n_mat = create_unique_mat();
  PetscErrorCode perr = MatDuplicate(matrix, MAT_DO_NOT_COPY_VALUES, e_n_mat.get());
  CHKERRXX(perr);
  perr = MatShift(*e_n_mat, -eigen_lo);
  CHKERRXX(perr);

  floating eigen_hi = get_eigenvalue_by_poweriter(*e_n_mat, conv_thres, max_iter) + eigen_lo;

  std::pair<floating, floating> eigenpair = std::minmax(std::abs(eigen_lo), std::abs(eigen_hi));

  eigen_lo = eigenpair.first;
  eigen_hi = eigenpair.second;
  floating condnum = eigen_hi/eigen_lo;

  fedisableexcept(FE_INVALID | FE_DIVBYZERO);

  return condnum; 
}

floating get_eigenvalue_by_poweriter(const Mat& matrix, floating conv_thres, integer max_iter)
{
  if (conv_thres <= 0)
  {
    throw InternalError("conv_thres must be greater than 0", __FILE__, __LINE__);
  }

  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);

  // Initial guess at b_0
  Vec_unique b = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(matrix, b.get(), nullptr);
  CHKERRXX(perr);
  perr = VecSet(*b, 1.0);
  CHKERRXX(perr);

  Vec_unique b_next = create_unique_vec();
  perr = MatCreateVecs(matrix, b_next.get(), nullptr);
  CHKERRXX(perr);

  integer curr_iter = 0;
  floating delta = 2 * conv_thres;
  floating lambda = 1.0;
  floating b_dot;
  perr = VecDot(*b, *b, &b_dot);
  CHKERRXX(perr);

  while (delta > conv_thres && curr_iter++ < max_iter)
  {
    perr = MatMult(matrix, *b, *b_next);
    CHKERRXX(perr);

    floating b_bn_dot;
    perr = VecDot(*b, *b_next, &b_bn_dot);
    CHKERRXX(perr);

    floating new_lambda = b_bn_dot / b_dot;
    delta = std::abs(new_lambda - lambda);

    perr = VecDot(*b_next, *b_next, &b_dot);
    CHKERRXX(perr);
    perr = VecScale(*b_next, 1. / std::sqrt(b_dot));
    CHKERRXX(perr);

    b.swap(b_next);
    lambda = new_lambda;
  }

  if (curr_iter >= max_iter)
  {
    PetscPrintf(comm, "Warning: Power iteration did not converge, "
                      "eigenvalue approximation may be inaccurate.");
  }

  return lambda;
}


void repeat_stack(const Vec& subvec, const Vec& stacked)
{
  MPI_Comm comm, comm2;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(subvec), &comm);
  PetscObjectGetComm(reinterpret_cast<PetscObject>(stacked), &comm2);
  if (comm != comm2){
    throw InternalError("stacked and subvec must share a communicator", __FILE__, __LINE__);
  }

  integer subvec_len;
  VecGetSize(subvec, &subvec_len);
  integer stacked_len;
  VecGetSize(stacked, &stacked_len);

  if (stacked_len % subvec_len != 0){
    std::ostringstream errss;
    errss << "stacked length (" << stacked_len << ") not a multiple of subvec length ("
          << subvec_len << ").";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  // Get extents of local data in subvector
  integer startelem, datasize;
  PetscErrorCode perr = VecGetOwnershipRange(subvec, &startelem, &datasize);
  CHKERRXX(perr);
  datasize -= startelem;
  // Create IS for source indices in stacked array
  IS_unique src_is = create_unique_is();
  perr = ISCreateStride(comm, datasize, startelem, 1, src_is.get());
  CHKERRXX(perr);

  integer stack_total = stacked_len / subvec_len;
  integer offset = 0; // Track offset
  for (integer count = 0; count < stack_total; count++)
  {
    
    IS_unique tgt_is = create_unique_is();
    perr = ISCreateStride(comm, datasize, startelem + offset, 1, tgt_is.get());
    CHKERRXX(perr);

    // Create scatterer and add to array
    VecScatter_unique scatterer = create_unique_vecscatter();
    perr = VecScatterCreate(subvec, *src_is, stacked, *tgt_is, scatterer.get());
    CHKERRXX(perr);

    perr = VecScatterBegin(*scatterer, subvec, stacked, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRXX(perr);
    perr = VecScatterEnd(*scatterer, subvec, stacked, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRXX(perr);

    offset += subvec_len;
  }
}


void repeat_stack_petsc_to_nat(const Vec& subvec, const Vec& stacked, const DM& subvec_dmda)
{
  MPI_Comm comm, comm2;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(subvec), &comm);
  PetscObjectGetComm(reinterpret_cast<PetscObject>(stacked), &comm2);
  if (comm != comm2){
    throw InternalError("stacked and subvec must share a communicator", __FILE__, __LINE__);
  }

  integer subvec_len;
  VecGetSize(subvec, &subvec_len);
  integer stacked_len;
  VecGetSize(stacked, &stacked_len);

  if (stacked_len % subvec_len != 0){
    std::ostringstream errss;
    errss << "stacked length (" << stacked_len << ") not a multiple of subvec length ("
          << subvec_len << ").";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed
                    // obj.
  PetscErrorCode perr = DMDAGetAO(subvec_dmda, &ao_petsctonat); // Destroying this would break the dmda
  CHKERRXX(perr);

  // Get extents of local data in subvector
  integer startelem, datasize;
  perr = VecGetOwnershipRange(subvec, &startelem, &datasize);
  CHKERRXX(perr);
  datasize -= startelem;
  // Create IS for source indices in stacked array
  IS_unique src_is = create_unique_is();
  perr = ISCreateStride(comm, datasize, startelem, 1, src_is.get());
  CHKERRXX(perr);
  // Convert with AO to map from natural to petsc ordering (do here because tgt_is is offset)
  perr = AOApplicationToPetscIS(ao_petsctonat, *src_is);
  CHKERRXX(perr);

  integer stack_total = stacked_len / subvec_len;
  integer offset = 0; // Track offset
  for (integer count = 0; count < stack_total; count++)
  {
    
    IS_unique tgt_is = create_unique_is();
    perr = ISCreateStride(comm, datasize, startelem + offset, 1, tgt_is.get());
    CHKERRXX(perr);

    // Create scatterer and add to array
    VecScatter_unique scatterer = create_unique_vecscatter();
    perr = VecScatterCreate(subvec, *src_is, stacked, *tgt_is, scatterer.get());
    CHKERRXX(perr);

    perr = VecScatterBegin(*scatterer, subvec, stacked, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRXX(perr);
    perr = VecScatterEnd(*scatterer, subvec, stacked, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRXX(perr);

    offset += subvec_len;
  }
}


void copy_nth_from_stack_nat_to_petsc(const Vec& subvec, const Vec& stacked, const DM& subvec_dmda, integer n)
{
  MPI_Comm comm, comm2;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(subvec), &comm);
  PetscObjectGetComm(reinterpret_cast<PetscObject>(stacked), &comm2);
  if (comm != comm2){
    throw InternalError("stacked and subvec must share a communicator", __FILE__, __LINE__);
  }

  integer subvec_len;
  VecGetSize(subvec, &subvec_len);
  integer stacked_len;
  VecGetSize(stacked, &stacked_len);
  integer stack_total = stacked_len / subvec_len;

  if (stacked_len % subvec_len != 0){
    std::ostringstream errss;
    errss << "stacked length (" << stacked_len << ") not a multiple of subvec length ("
          << subvec_len << ").";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  if (n >= stack_total)
  {
    std::ostringstream errss;
    errss << "n must be less than stacked length / subvec length.";
    throw InternalError(errss.str(), __FILE__, __LINE__);
  }

  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed
                    // obj.
  PetscErrorCode perr = DMDAGetAO(subvec_dmda, &ao_petsctonat); // Destroying this would break the dmda
  CHKERRXX(perr);

  // Get extents of local data in subvector
  integer startelem, datasize;
  perr = VecGetOwnershipRange(subvec, &startelem, &datasize);
  CHKERRXX(perr);
  datasize -= startelem;
  // Create IS for source indices in stacked array
  IS_unique src_is = create_unique_is();
  perr = ISCreateStride(comm, datasize, startelem, 1, src_is.get());
  CHKERRXX(perr);
  // Convert with AO to map from natural to petsc ordering (do here because tgt_is is offset)
  // This direction because then we'll reverse the scatter
  perr = AOApplicationToPetscIS(ao_petsctonat, *src_is);
  CHKERRXX(perr);

  integer offset = n*subvec_len; // Track offset
  IS_unique tgt_is = create_unique_is();
  perr = ISCreateStride(comm, datasize, startelem + offset, 1, tgt_is.get());
  CHKERRXX(perr);

  // Create scatterer and add to array
  VecScatter_unique scatterer = create_unique_vecscatter();
  perr = VecScatterCreate(subvec, *src_is, stacked, *tgt_is, scatterer.get());
  CHKERRXX(perr);

  perr = VecScatterBegin(*scatterer, stacked, subvec, INSERT_VALUES, SCATTER_REVERSE);
  CHKERRXX(perr);
  perr = VecScatterEnd(*scatterer, stacked, subvec, INSERT_VALUES, SCATTER_REVERSE);
  CHKERRXX(perr);
}


