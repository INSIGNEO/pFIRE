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

void mainflow();

int main(int argc, char **argv){
  PetscErrorCode perr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRABORT(PETSC_COMM_WORLD, perr);

  mainflow();

  PetscFinalize();

  return 0;
}

void mainflow(){

  intvector imshape = {3,3,3};
  intvector nodespacing = {2,2,2};

  Image testimg(imshape);

  Map testmap(nodespacing, testimg);

  MatView(*(testmap.basis()), PETSC_VIEWER_STDOUT_WORLD);

}
