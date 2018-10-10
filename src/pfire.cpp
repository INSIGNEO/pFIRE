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

  floatvector nodespacing = {5,5,5};

  std::unique_ptr<Image> fixed = Image::create_from_image("ellipse_fixed.png");
  std::unique_ptr<Image> moved = Image::create_from_image("ellipse_moved.png", fixed.get());

  Elastic reg(*fixed, *moved, nodespacing);
  reg.autoregister();
}
