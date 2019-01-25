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

bool vecs_equivalent(const Vec &vec1, const Vec &vec2)
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
