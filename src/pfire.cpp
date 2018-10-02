#include<petscmat.h>
#include<petscvec.h>
#include<petscdmda.h>
#include<petscviewer.h>
#include<mpi.h>

#include<iostream>
#include "types.hpp"
#include "laplacian.hpp"
#include "image.hpp"
#include "map.hpp"
#include "elastic.hpp"

void mainflow();

int main(int argc, char **argv){
  PetscErrorCode perr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRABORT(PETSC_COMM_WORLD, perr);

  mainflow();

  PetscFinalize();

  return 0;
}

void mainflow(){

  intvector imshape = {10,10,10};
  floatvector nodespacing = {6,6,6};
  floatvector nodespacing2 = {3,3,3};

  std::unique_ptr<Image> fixed = std::make_unique<Image>(imshape);
  std::unique_ptr<Image> moved = fixed->duplicate();

  Map foo(*fixed, nodespacing);

  VecSet(*foo.m_displacements, 2);

  auto foo2 = foo.interpolate(nodespacing2);

  VecView(*foo2->m_displacements, PETSC_VIEWER_STDOUT_WORLD);

}
