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

  intvector imshape = {10,10};
  floatvector nodespacing = {5,5};

  std::unique_ptr<Image> fixed = std::make_unique<Image>(imshape);
  std::unique_ptr<Image> moved = fixed->duplicate();

  Elastic reg(*fixed, *moved, nodespacing);
  reg.autoregister();
}
