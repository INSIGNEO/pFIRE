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

floating diagonal_sum(const Mat& matrix)
{
  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);

  Vec_unique diag = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(matrix, diag.get(), nullptr);
  CHKERRABORT(comm, perr);

  perr = MatGetDiagonal(matrix, *diag);
  CHKERRABORT(comm, perr);

  floating diagsum;
  perr = VecSum(*diag, &diagsum);

  return diagsum;
}

bool vecs_equivalent(const Vec& vec1, const Vec& vec2)
{
  // Ensure two vectors have the same local size, global size. and comm.
  MPI_Comm comm1, comm2;
  PetscErrorCode perr = PetscObjectGetComm((PetscObject)vec1, &comm1);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = PetscObjectGetComm((PetscObject)vec2, &comm2);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (comm1 != comm2)
  {
    return false;
  }

  integer globalsize1, globalsize2;
  perr = VecGetSize(vec1, &globalsize1);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecGetSize(vec2, &globalsize2);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (globalsize1 != globalsize2)
  {
    return false;
  }

  integer localsize1, localsize2;
  perr = VecGetLocalSize(vec1, &localsize1);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecGetLocalSize(vec2, &localsize2);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
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
  CHKERRABORT(comm, perr);
  debug_creation(*diag, "Vec_block_normalise");
  perr = MatGetDiagonal(normmat, *diag);
  CHKERRABORT(comm, perr);

  // index where spatial dims stop and intensity dim starts
  integer crit_idx = size * ndim;

  // Find rank-local sums for spatial and luminance blocks
  integer rowstart, rowend;
  perr = VecGetOwnershipRange(*diag, &rowstart, &rowend);
  CHKERRABORT(comm, perr);
  integer localsize = rowend - rowstart;
  floatvector norm(2, 0.0); // norm[0] is spatial, norm[1] is luminance
  integer spt_end = crit_idx - rowstart;
  spt_end = (spt_end < 0) ? 0 : (spt_end > localsize) ? localsize : spt_end;

  floating* ptr;
  perr = VecGetArray(*diag, &ptr);
  CHKERRABORT(comm, perr);
  for (integer idx = 0; idx < spt_end; idx++)
  {
    norm[0] += ptr[idx];
  }
  for (integer idx = spt_end; idx < localsize; idx++)
  {
    norm[1] += ptr[idx];
  }
  perr = VecRestoreArray(*diag, &ptr);
  CHKERRABORT(comm, perr);

  // MPI_AllReduce to sum over all processes
  MPI_Allreduce(MPI_IN_PLACE, norm.data(), 2, MPI_DOUBLE, MPI_SUM, comm);

  // calculate average of norms and scaling factor
  norm[0] /= crit_idx;
  norm[1] /= size;
  floating lum_scale = norm[0] / norm[1];

  // reuse diag vector to hold scaling values
  perr = VecGetArray(*diag, &ptr);
  CHKERRABORT(comm, perr);
  for (integer idx = 0; idx < spt_end; idx++)
  {
    ptr[idx] = 1;
  }
  for (integer idx = spt_end; idx < localsize; idx++)
  {
    ptr[idx] = lum_scale;
  }
  perr = VecRestoreArray(*diag, &ptr);
  CHKERRABORT(comm, perr);

  perr = MatDiagonalScale(normmat, *diag, nullptr);
  CHKERRABORT(comm, perr);
}

floating get_condnum_by_poweriter(const Mat& matrix, floating conv_thres, integer max_iter)
{
  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);
  floating eigen_lo = get_eigenvalue_by_poweriter(matrix, conv_thres, max_iter);

  // Get smallest eigenvalue by spectral shift
  Mat_unique e_n_mat = create_unique_mat();
  PetscErrorCode perr = MatDuplicate(matrix, MAT_DO_NOT_COPY_VALUES, e_n_mat.get());
  CHKERRABORT(comm, perr);
  perr = MatShift(*e_n_mat, -eigen_lo);
  CHKERRABORT(comm, perr);

  floating eigen_hi = get_eigenvalue_by_poweriter(*e_n_mat, conv_thres, max_iter) + eigen_lo;

  std::pair<floating, floating> eigenpair = std::minmax(std::abs(eigen_lo), std::abs(eigen_hi));

  eigen_lo = eigenpair.first;
  eigen_hi = eigenpair.second;

  return eigen_hi / eigen_lo;
}

floating get_eigenvalue_by_poweriter(const Mat& matrix, floating conv_thres, integer max_iter)
{
  if (conv_thres <= 0)
  {
    throw std::domain_error("conv_thres must be greater than 0");
  }

  MPI_Comm comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);

  // Initial guess at b_0
  Vec_unique b = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(matrix, b.get(), nullptr);
  CHKERRABORT(comm, perr);
  perr = VecSet(*b, 1.0);
  CHKERRABORT(comm, perr);

  Vec_unique b_next = create_unique_vec();
  perr = MatCreateVecs(matrix, b_next.get(), nullptr);
  CHKERRABORT(comm, perr);

  integer curr_iter = 0;
  floating delta = 2 * conv_thres;
  floating lambda = 1.0;
  floating b_dot;
  perr = VecDot(*b, *b, &b_dot);
  CHKERRABORT(comm, perr);

  while (delta > conv_thres && curr_iter++ < max_iter)
  {
    perr = MatMult(matrix, *b, *b_next);
    CHKERRABORT(comm, perr);

    floating b_bn_dot;
    perr = VecDot(*b, *b_next, &b_bn_dot);
    CHKERRABORT(comm, perr);

    floating new_lambda = b_bn_dot / b_dot;
    delta = std::abs(new_lambda - lambda);

    perr = VecDot(*b_next, *b_next, &b_dot);
    CHKERRABORT(comm, perr);
    perr = VecScale(*b_next, 1. / std::sqrt(b_dot));
    CHKERRABORT(comm, perr);

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

// Determine suitable grid decomposition by minimizing chunk surface area
std::pair<intvector, intvector> find_proc_split(
    const intvector& griddims, const integer& comm_size)
{
  intvector gridsplit(3, 0);
  intvector procsplit(3, 0);

  integer minarea =
      griddims[0] * griddims[1] + griddims[1] * griddims[2] + griddims[2] * griddims[0];

  integer zcomm_eff_size = comm_size;
  if (griddims[2] == 1)
  {
    zcomm_eff_size = 1;
  }
  // Iterate over the comm and find suitable size
  for (integer iz = 1; iz <= zcomm_eff_size; iz++)
  {
    integer nprocxy = comm_size / iz;
    if (iz * nprocxy != comm_size)
      continue;

    // If the splitting is too fine no point in increasing rank in z dir
    integer zsplit = griddims[2] / iz;
    if (zsplit < 1 && iz > 1)
      break;

    // Now find y and x size
    for (integer iy = 1; iy <= nprocxy; iy++)
    {
      integer ix = nprocxy / iy;
      if (iy * ix != nprocxy)
        continue;

      integer ysplit = griddims[1] / iy;
      integer xsplit = griddims[0] / ix;

      // Skip if splitting too fine
      if (xsplit < 1 || ysplit < 1)
        continue;

      // If best yet, save result as potential winner
      integer area = xsplit * ysplit + ysplit * zsplit + zsplit * xsplit;
      if (area <= minarea)
      {
        gridsplit[0] = xsplit;
        gridsplit[1] = ysplit;
        gridsplit[2] = zsplit;
        procsplit[0] = ix;
        procsplit[1] = iy;
        procsplit[2] = iz;
      }
    }
  }

  if (gridsplit[0] > 0)
  {
    return std::make_pair(procsplit, gridsplit);
  }

  // We could *try* making a new comm here, but it makes a mess of everything else
  throw std::runtime_error("Unable to partition image appropriately");
}
