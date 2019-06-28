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

#include "tmatrix.hpp"

#include <cfenv>
#include <chrono>
#include <iostream>
#include <iterator>
#include <thread>

#include "fd_routines.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "infix_iterator.hpp"
#include "mapbase.hpp"
#include "mask.hpp"
#include "mpi_routines.hpp"

/*! Build the T^t T matrix object and T^f(f-m) vector object
 *
 */
std::pair<Mat_unique, Vec_unique> build_tmat2_and_tmatfm(const Image& fixed, const Image& moved,
                                                         const MapBase& map)
{
  MPI_Comm comm = fixed.comm();
  int commsize, rank;
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &rank);

  // Final matrix dimensions follow from the map
  integer global_rows = map.num_total_nodes() * map.ndof();

  // Calculate gradient data for later
  auto gradient_data = calculate_tmatrix_gradients(fixed, moved, map.ndim());

  Vec_unique image_difference = fixed.difference_global(moved);

  // Calculate offsets for columns that produce non-zero entries with given row
  auto col_offsets = calculate_moore_offsets(map.shape(), map.ndof());

  // Create T^t(f-m) vector storage
  Vec_unique tt_fm = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(map.displacement_vector(), tt_fm.get());
  CHKERRXX(perr);
  perr = VecSet(*tt_fm, 0.);
  CHKERRXX(perr);

  // Find our row range matching the tt_fm entries
  integer row_lo, row_hi;
  perr = VecGetOwnershipRange(*tt_fm, &row_lo, &row_hi);
  CHKERRXX(perr);
  integer local_rows = row_hi - row_lo;

  // Reserve data storage for array elements
  intvector row_pointers;
  row_pointers.reserve(local_rows + 1);
  intvector column_pointers;
  column_pointers.reserve(local_rows * map.stencilpoints()); // Can have up to 21 entries per row
  floatvector matrix_data;
  matrix_data.reserve(local_rows * map.stencilpoints());
  floatvector vec_entries;
  vec_entries.reserve(local_rows * map.ndof());

  // communicate max interations
  integer max_iter = local_rows + map.ndof();
  MPI_Allreduce(MPI_IN_PLACE, &max_iter, 1, MPIU_INT, MPI_MAX, comm);

  // Iterate over all our owned rows
  integer iter_count = 0;
  integer data_pointer(0);
  integer dummy_row_lo = row_lo - (row_lo % map.ndof());

  // Start from idof=0 to keep all ranks working on same data
  for (integer curr_row = dummy_row_lo; curr_row < row_lo; curr_row++)
  {
    integer idof = curr_row % map.ndof();
    for (const auto& col_offset __attribute__((unused)) : col_offsets)
    {
      communicate_entries_only(idof, map, gradient_data, *image_difference, comm);
    }
    iter_count++;
  }

  for (integer curr_row = row_lo; curr_row < row_hi; curr_row++)
  {
    floating vec_entry(0);
    integer idof = curr_row % map.ndof();
    row_pointers.push_back(data_pointer);
    // Iterate over possible non-zero columns for mult
    for (const auto& col_offset : col_offsets)
    {
      integer curr_col = curr_row + col_offset; // account for ndof in map
      if (curr_col < 0 || curr_col >= global_rows)
      {
        communicate_entries_only(idof, map, gradient_data, *image_difference, comm);
        continue;
      }
      auto entries = calculate_entries(curr_row, curr_col, idof, map, gradient_data, *image_difference,
                                       comm);
      vec_entry += entries.second;

      floating& mat_entry = entries.first;
      column_pointers.push_back(curr_col);
      matrix_data.push_back(mat_entry);
      data_pointer++;
    }
    vec_entries.push_back(vec_entry);
    iter_count++;
  }
  row_pointers.push_back(data_pointer);

  while (iter_count < max_iter)
  {
    integer idof = iter_count % map.ndof();
    for (const auto& col_offset __attribute__((unused)) : col_offsets)
    {
      communicate_entries_only(idof, map, gradient_data, *image_difference, comm);
    }
    iter_count++;
  }

  // Insert vector entries
  floating* vecarr;
  perr = VecGetArray(*tt_fm, &vecarr);
  CHKERRXX(perr);
  std::copy(vec_entries.cbegin(), vec_entries.cend(), vecarr);
  perr = VecRestoreArray(*tt_fm, &vecarr);
  CHKERRXX(perr);

  // Finally construct our CSR matrix
  Mat_unique tmat2 = create_unique_mat();
  perr = MatCreateMPIAIJWithArrays(comm, local_rows,
                                   local_rows,               // Local matrix size
                                   global_rows, global_rows, // Specify as cross-check
                                   row_pointers.data(),      // CSR data
                                   column_pointers.data(),   //   "
                                   matrix_data.data(),       //   "
                                   tmat2.get());
  CHKERRXX(perr);

  precondition_tmat2(*tmat2, map.ndof(), comm);

  return std::make_pair(std::move(tmat2), std::move(tt_fm));
}

void communicate_entries_only(integer idof, const MapBase& map, const std::vector<Vec_unique>& gradient_data,
                              const Vec& difference_vec, MPI_Comm comm)
{
#ifdef DEBUG_CHECKS
  throw_if_idof_inconsistent(comm, idof);
#endif //DEBUG_CHECKS

  int commsize;
  MPI_Comm_size(comm, &commsize);

#ifdef VERBOSE_DEBUG
  for (int irank = 0; irank < commsize; irank++)
  {
    MPI_Barrier(comm);
  }
#endif // VERBOSE_DEBUG

  intcoordvector2d remote_reqs(commsize);
  // scatter as needed
  if (commsize > 1)
  {
    remote_reqs = p2p_vecscatter(remote_reqs, comm);
  }

  // look up the values
  floatpairvector2d remote_vals = populate_request_vector(remote_reqs, map.mask().dmda(),
                                                          *gradient_data[idof], difference_vec);

  // scatter back as needed
  if (commsize > 1)
  {
    remote_vals = p2p_vecscatter(remote_vals, comm);
  }

#ifdef VERBOSE_DEBUG
  for (int irank = 0; irank < commsize; irank++)
  {
    MPI_Barrier(comm);
  }

  for (int irank = 0; irank < commsize; irank++)
  {
    MPI_Barrier(comm);
  }
#endif // VERBOSE_DEBUG
}

std::pair<floating, floating> calculate_entries(integer row, integer col, integer idof, const MapBase& map,
                                                const std::vector<Vec_unique>& gradient_data,
                                                const Vec& difference_data, MPI_Comm comm)
{
#ifdef DEBUG_CHECKS
  throw_if_idof_inconsistent(comm, idof)
#endif //DEBUG_CHECKS

  int commsize, rank;
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &rank);

#ifdef DEBUG_CHECKS
  if (idof != row % map.ndof() || idof != col % map.ndof())
  {
    throw InternalError("Wrong dof in calculation");
  }
#endif //DEBUG_CHECKS

  // First calculate appropriate overlaps
  auto rloc = unravel(static_cast<integer>(row / map.ndof()), map.shape());
  auto cloc = unravel(static_cast<integer>(col / map.ndof()), map.shape());

  auto row_pixrange = map.get_pixel_neighbourhood(rloc); // returns pairs of (lo, hi) intcoords
  auto col_pixrange = map.get_pixel_neighbourhood(cloc);

  auto overlap_range = get_overlap(row_pixrange.first, row_pixrange.second, // another range pair
                                   col_pixrange.first, col_pixrange.second);

  intcoord overlap_curr = overlap_range.first;
  auto& xx = overlap_curr[0];
  auto& yy = overlap_curr[1];
  auto& zz = overlap_curr[2];

  // Create comm temporary storage
  floatvector2d row_coeffs(commsize);
  floatvector2d col_coeffs(commsize);
  intcoordvector2d remote_reqs(commsize);

#ifdef VERBOSE_DEBUG
  for (int irank = 0; irank < commsize; irank++)
  {
    if (rank == irank)
    {
      std::cout << "(" << row << ", " << col << ") -> " << overlap_range.first << " -- "
                << overlap_range.second << "\n"
                << std::flush;
    }
    MPI_Barrier(comm);
  }
#endif // VERBOSE_DEBUG

  for (xx = overlap_range.first[0]; xx < overlap_range.second[0]; xx++)
  {
    for (yy = overlap_range.first[1]; yy < overlap_range.second[1]; yy++)
    {
      for (zz = overlap_range.first[2]; zz < overlap_range.second[2]; zz++)
      {
        floating r_coeff = calculate_scaled_basis_coefficient(overlap_curr.cbegin(), overlap_curr.cend(),
                                                              map.coord_from_index(rloc).cbegin(),
                                                              map.spacing().cbegin());
        floating c_coeff = calculate_scaled_basis_coefficient(overlap_curr.cbegin(), overlap_curr.cend(),
                                                              map.coord_from_index(cloc).cbegin(),
                                                              map.spacing().cbegin());
        // lookup remote rank
        integer req_rank = map.mask().get_rank_of_loc(overlap_curr);
        // insert row/column_coeffs
        row_coeffs[req_rank].push_back(r_coeff);
        col_coeffs[req_rank].push_back(c_coeff);
        // add to request structure
        remote_reqs[req_rank].push_back(overlap_curr);
      }
    }
  }

  // scatter as needed
  //  if (commsize > 1)
  //  {
  auto commed_remote_reqs = p2p_vecscatter(remote_reqs, comm);
  //  }

  // look up the values
  floatpairvector2d uncommed_remote_vals = populate_request_vector(commed_remote_reqs, map.mask().dmda(),
                                                                   *gradient_data[idof], difference_data);

  // scatter back as needed
  //  if (commsize > 1)
  //  {
  auto remote_vals = p2p_vecscatter(uncommed_remote_vals, comm);
  //  }

  floating mat_entry(0);
  floating mat_err(0);
  floating vec_entry(0);
  floating vec_err(0);
#ifdef VERBOSE_DEBUG
  integer sumcount(0);
  for (int irank = 0; irank < commsize; irank++)
  {
    if (irank == rank)
    {
#endif // VERBOSE_DEBUG
      for (int jrank = 0; jrank < commsize; jrank++)
      {
#ifdef VERBOSE_DEBUG
        auto& rank_reqs = remote_reqs[jrank];
#endif // VERBOSE_DEBUG
        auto& rank_vals = remote_vals[jrank];
        auto& rank_row_coeffs = row_coeffs[jrank];
        auto& rank_col_coeffs = col_coeffs[jrank];
        for (size_t iidx = 0; iidx < rank_vals.size(); iidx++)
        {
          floating mat_delta = (rank_vals[iidx].first * rank_row_coeffs[iidx])
                               * (rank_vals[iidx].first * rank_col_coeffs[iidx]);
          mat_delta -= mat_err;
          floating tmp = mat_entry + mat_delta;
          mat_err = (tmp - mat_entry) - mat_delta;
          mat_entry = tmp;

          floating vec_delta = 0.5 * rank_vals[iidx].first
                               * (rank_vals[iidx].second * rank_row_coeffs[iidx]);
          vec_delta -= vec_err;
          tmp = vec_entry + vec_delta;
          vec_err = (tmp - vec_entry) - vec_delta;
          vec_entry = tmp;

#ifdef VERBOSE_DEBUG
          sumcount++;
          std::cout << "Component (" << row << ", " << col << ") " << rank_reqs[iidx] << "[" << idof << "]"
                    << ", row_coeff = " << rank_row_coeffs[iidx] << ", col_coeff = " << rank_col_coeffs[iidx]
                    << ", grad = " << rank_vals[iidx].first << ", diff = " << rank_vals[iidx].second
                    << ", mat_delta = " << mat_delta << ", vec_delta = " << vec_delta << std::endl;
#endif // VERBOSE_DEBUG
        }
      }
#ifdef VERBOSE_DEBUG
      std::cout << std::flush;
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    MPI_Barrier(comm);
  }
#endif // VERBOSE_DEBUG

#ifdef VERBOSE_DEBUG
  for (int irank = 0; irank < commsize; irank++)
  {
    if (irank == rank)
    {
      std::cout << "Entry (" << row << ", " << col << ") sc: " << sumcount << ", me:" << mat_entry
                << " ve: " << vec_entry << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    MPI_Barrier(comm);
  }
#endif // VERBOSE_DEBUG

  return std::make_pair(mat_entry, vec_entry);
}

/*! Return remotely requested data from local vector
 *
 */
floatpairvector2d populate_request_vector(intcoordvector2d remote_locs, const DM& dmda,
                                          const Vec& gradient_vec, const Vec& difference_vec)
{
  floatpairvector2d remote_vals(remote_locs.size());

  floating const ***gradient_data, ***difference_data;
  PetscErrorCode perr = DMDAVecGetArrayRead(dmda, gradient_vec, &gradient_data);
  CHKERRXX(perr);
  perr = DMDAVecGetArrayRead(dmda, difference_vec, &difference_data);
  CHKERRXX(perr);

  DMDALocalInfo dminfo;
  DMDAGetLocalInfo(dmda, &dminfo);
  intcoord lo = {dminfo.xs, dminfo.ys, dminfo.zs};
  intcoord hi = {dminfo.xm, dminfo.ym, dminfo.zm};
  hi = hi + lo;
  

  for (size_t irank = 0; irank < remote_locs.size(); irank++)
  {
    auto& icoeffs = remote_vals[irank];
    auto& ilocs = remote_locs[irank];
    icoeffs.resize(ilocs.size());
    for (size_t idx = 0; idx < ilocs.size(); idx++)
    {
      auto& iloc = ilocs[idx];
#ifdef DEBUG_CHECKS
      if(iloc < lo || iloc >= hi)
      {
        throw InternalError("Loc not on rank");
      }
#endif //DEBUG_CHECKS
      icoeffs[idx] = std::make_pair(gradient_data[iloc[2]][iloc[1]][iloc[0]],
                                    difference_data[iloc[2]][iloc[1]][iloc[0]]);
    }
  }
  perr = DMDAVecRestoreArrayRead(dmda, gradient_vec, &gradient_data);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArrayRead(dmda, difference_vec, &difference_data);
  CHKERRXX(perr);

  return remote_vals;
}

void precondition_tmat2(const Mat& tmat2, integer ndof, MPI_Comm _comm)
{
  Vec_unique diag;
  diag = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(tmat2, diag.get(), nullptr);
  CHKERRXX(perr);
  perr = MatGetDiagonal(tmat2, *diag);
  CHKERRXX(perr);

  integer vec_size;
  perr = VecGetSize(*diag, &vec_size);

  integer row_lo, row_hi;
  perr = VecGetOwnershipRange(*diag, &row_lo, &row_hi);
  CHKERRXX(perr);

  std::vector<floating> sumdata(2);
  floating& sum1 = sumdata[0];
  floating& sum2 = sumdata[1];
  floating* diag_data;
  perr = VecGetArray(*diag, &diag_data);
  CHKERRXX(perr);
  for (integer iidx = 0; iidx < (row_hi - row_lo); iidx++)
  {
    if ((iidx + row_lo) % ndof == 0)
    {
      sum1 += diag_data[iidx];
      continue;
    }
    sum2 += diag_data[iidx];
  }

  MPI_Allreduce(MPI_IN_PLACE, sumdata.data(), sumdata.size(), MPIU_SCALAR, MPI_SUM, _comm);

  sum1 /= double(vec_size) / ndof;
  sum2 /= (ndof - 1) * double(vec_size) / ndof;
  floating mult = sum2 / sum1;

  for (integer iidx = 0; iidx < (row_hi - row_lo); iidx++)
  {
    if ((iidx + row_lo) % ndof == 0)
    {
      diag_data[iidx] = mult;
      continue;
    }
    diag_data[iidx] = 1;
  }
  perr = VecRestoreArray(*diag, &diag_data);
  CHKERRXX(perr);

  perr = MatDiagonalScale(tmat2, *diag, nullptr);
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
  floating condnum = eigen_hi / eigen_lo;

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

void throw_if_idof_inconsistent(MPI_Comm comm, integer idof)
{
  int commsize;
  MPI_Comm_size(comm, &commsize);

  std::vector<integer> idofs(commsize, idof);

  MPI_Alltoall(MPI_IN_PLACE, 1, MPIU_INT, idofs.data(), 1, MPIU_INT, comm);

  if(std::accumulate(idofs.begin(), idofs.end(), 0) != (idof*commsize))
  {
    std::ostringstream errss;
    errss << "idofs inconsistent: ";
    std::copy(idofs.begin(), idofs.end(), infix_ostream_iterator<integer>(errss, " "));
    throw InternalError(errss.str());
  }

}
