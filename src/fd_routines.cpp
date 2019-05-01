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

Vec_unique fd::gradient_to_global_unique(const DM &dmda, const Vec &localvec, integer dim)
{
  //  First sanity check we have a valid local vector for the DMDA
  //  (should have matching global, local and comm size)
  Vec dm_local_vec;
  PetscErrorCode perr = DMGetLocalVector(dmda, &dm_local_vec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (!vecs_equivalent(dm_local_vec, localvec))
  {
    throw InternalError("provided vector invalid for given dmda object", __FILE__, __LINE__);
  }
  perr = DMRestoreLocalVector(dmda, &dm_local_vec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  // New global vec must be a dmda global
  Vec_unique grad = create_unique_vec();
  perr = DMCreateGlobalVector(dmda, grad.get());
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  // Access local of this to have ghosts, global of grad to avoid later copy
  // Recall grad and image share a DMDA
  floating ***img_array,
      ***grad_array; // acceptable to use raw ptrs here as memory handled by PETSc
  perr = DMDAVecGetArray(dmda, localvec, &img_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMDAVecGetArray(dmda, *(grad.get()), &grad_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  integer i_lo, i_hi, j_lo, j_hi, k_lo, k_hi;
  // This returns corners + widths so add lo to get hi
  perr = DMDAGetCorners(dmda, &i_lo, &j_lo, &k_lo, &i_hi, &j_hi, &k_hi);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
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
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMDAVecRestoreArray(dmda, *(grad.get()), &grad_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  return grad;
}

void fd::gradient_existing(const DM &dmda, const Vec &srcvec, Vec &tgtvec, integer dim)
{
  //  First sanity check we have valid local vectors for the DMDA
  //  (should have matching global, local and comm size)
  Vec dm_local_vec;
  PetscErrorCode perr = DMGetLocalVector(dmda, &dm_local_vec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (!vecs_equivalent(dm_local_vec, srcvec))
  {
    throw InternalError("provided srcvec invalid for given dmda object", __FILE__, __LINE__);
  }
  perr = DMRestoreLocalVector(dmda, &dm_local_vec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  // Access local of this to have ghosts, global of grad to avoid later copy
  // Recall grad and image share a DMDA
  floating ***img_array,
      ***grad_array; // acceptable to use raw ptrs here as memory handled by PETSc
  perr = DMDAVecGetArray(dmda, srcvec, &img_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMDAVecGetArray(dmda, tgtvec, &grad_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  integer i_lo, i_hi, j_lo, j_hi, k_lo, k_hi;
  // This returns corners + widths so add lo to get hi
  perr = DMDAGetCorners(dmda, &i_lo, &j_lo, &k_lo, &i_hi, &j_hi, &k_hi);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
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
        // Remember c-indexing is backwards because PETSc is odd
        grad_array[k][j][i] = 0.5
                              * (img_array[k + ofs[2]][j + ofs[1]][i + ofs[0]]
                                 - img_array[k - ofs[2]][j - ofs[1]][i - ofs[0]]);
      }
    }
  }
  // Release pointers and allow petsc to cleanup
  perr = DMDAVecRestoreArray(dmda, srcvec, &img_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMDAVecRestoreArray(dmda, tgtvec, &grad_array);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
}
