#include "petsc_helpers.hpp"

bool vecs_equivalent(const Vec &vec1, const Vec &vec2){
  //Ensure two vectors have the same local size, global size. and comm.
  MPI_Comm comm1, comm2;
  PetscErrorCode perr = PetscObjectGetComm((PetscObject)vec1, &comm1);CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = PetscObjectGetComm((PetscObject)vec2, &comm2);CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (comm1 != comm2) {return false;}
  
  integer globalsize1, globalsize2;
  perr = VecGetSize(vec1, &globalsize1);CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecGetSize(vec2, &globalsize2);CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (globalsize1 != globalsize2) {return false;}

  integer localsize1, localsize2;
  perr = VecGetLocalSize(vec1, &localsize1);CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecGetLocalSize(vec2, &localsize2);CHKERRABORT(PETSC_COMM_WORLD, perr);
  if (localsize1 != localsize2) {return false;}

  return true;
}


