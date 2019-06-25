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

#include "fd_routines.hpp"

#include "exceptions.hpp"

Vec_unique gradient_to_global_unique(const DM& dmda, const Vec& localvec, integer dim)
{
  // New global vec must be a dmda global
  Vec_unique grad = create_unique_vec();
  PetscErrorCode perr = DMCreateGlobalVector(dmda, grad.get());
  CHKERRXX(perr);

  // Access local of this to have ghosts, global of grad to avoid later copy
  // Recall grad and image share a DMDA
  floating ***img_array,
      ***grad_array; // acceptable to use raw ptrs here as memory handled by PETSc
  perr = DMDAVecGetArray(dmda, localvec, &img_array);
  CHKERRXX(perr);
  perr = DMDAVecGetArray(dmda, *(grad.get()), &grad_array);
  CHKERRXX(perr);

  integer i_lo, i_hi, j_lo, j_hi, k_lo, k_hi;
  // This returns corners + widths so add lo to get hi
  perr = DMDAGetCorners(dmda, &i_lo, &j_lo, &k_lo, &i_hi, &j_hi, &k_hi);
  CHKERRXX(perr);
  i_hi += i_lo;
  j_hi += j_lo;
  k_hi += k_lo;

  intvector ofs = {0, 0, 0};
  ofs[dim] = 1;
  for (integer i = i_lo; i < i_hi; i++)
  {
    for (integer j = j_lo; j < j_hi; j++)
    {
      for (integer k = k_lo; k < k_hi; k++)
      {
        // Remember c-indexing is backwards because PETSc is a fortran monster
        grad_array[k][j][i] = 0.5
                              * (img_array[k + ofs[2]][j + ofs[1]][i + ofs[0]]
                                 - img_array[k - ofs[2]][j - ofs[1]][i - ofs[0]]);
      }
    }
  }
  // Release pointers and allow petsc to cleanup
  perr = DMDAVecRestoreArray(dmda, localvec, &img_array);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArray(dmda, *(grad.get()), &grad_array);
  CHKERRXX(perr);

  return grad;
}

void gradient_existing(const DM& dmda, const Vec& srcvec, Vec& tgtvec, integer dim)
{
  // Access local of this to have ghosts, global of grad to avoid later copy
  // Recall grad and image share a DMDA
  floating ***img_array,
      ***grad_array; // acceptable to use raw ptrs here as memory handled by PETSc
  PetscErrorCode perr = DMDAVecGetArray(dmda, srcvec, &img_array);
  CHKERRXX(perr);
  perr = DMDAVecGetArray(dmda, tgtvec, &grad_array);
  CHKERRXX(perr);

  integer i_lo, i_hi, j_lo, j_hi, k_lo, k_hi;
  // This returns corners + widths so add lo to get hi
  perr = DMDAGetCorners(dmda, &i_lo, &j_lo, &k_lo, &i_hi, &j_hi, &k_hi);
  CHKERRXX(perr);
  i_hi += i_lo;
  j_hi += j_lo;
  k_hi += k_lo;

#ifdef GRADIENT_DEBUG
  MPI_Comm comm;
  int rank, commsize;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(dmda), &comm);
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &rank);
  for (int irank = 0; irank < commsize; irank++)
  {
    if (irank == rank)
    {
#endif // GRADIENT_DEBUG
      intcoord ofs = {0, 0, 0};
      ofs[dim] = 1;
      for (integer i = i_lo; i < i_hi; i++)
      {
        for (integer j = j_lo; j < j_hi; j++)
        {
          for (integer k = k_lo; k < k_hi; k++)
          {
#ifdef GRADIENT_DEBUG
            // Remember c-indexing is backwards because PETSc is odd
            std::cout << "(" << i << ", " << j << ", " << k << ") +/- " << ofs << ": "
                      << img_array[k + ofs[2]][j + ofs[1]][i + ofs[0]] << " - "
                      << img_array[k - ofs[2]][j - ofs[1]][i - ofs[0]] << std::endl;
#endif // GRADIENT_DEBUG
            grad_array[k][j][i] = 0.5
                                  * (img_array[k + ofs[2]][j + ofs[1]][i + ofs[0]]
                                     - img_array[k - ofs[2]][j - ofs[1]][i - ofs[0]]);
          }
        }
      }
#ifdef GRADIENT_DEBUG
    }
    MPI_Barrier(comm);
  }
#endif // GRADIENT_DEBUG

  // Release pointers and allow petsc to cleanup
  perr = DMDAVecRestoreArray(dmda, srcvec, &img_array);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArray(dmda, tgtvec, &grad_array);
  CHKERRXX(perr);

  perr = DMLocalToLocalBegin(dmda, tgtvec, INSERT_VALUES, tgtvec);
  CHKERRXX(perr);
  perr = DMLocalToLocalEnd(dmda, tgtvec, INSERT_VALUES, tgtvec);
  CHKERRXX(perr);
}
