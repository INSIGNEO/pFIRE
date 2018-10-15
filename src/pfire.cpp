#include<petscmat.h>
#include<petscvec.h>
#include<petscdmda.h>
#include<petscviewer.h>
#include<mpi.h>

#include<iostream>
#include<string>

#include "types.hpp"
#include "laplacian.hpp"
#include "image.hpp"
#include "map.hpp"
#include "elastic.hpp"

void mainflow(std::string, std::string);

int main(int argc, char **argv){
  PetscErrorCode perr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRABORT(PETSC_COMM_WORLD, perr);

  mainflow(argv[1], argv[2]);

  PetscFinalize();

  return 0;
}

void mainflow(std::string fixedpath, std::string movedpath){

  floatvector nodespacing = {5,5};

  std::unique_ptr<Image> fixed = Image::create_from_image(fixedpath);
  std::unique_ptr<Image> moved = Image::create_from_image(movedpath, fixed.get());

  Elastic reg(*fixed, *moved, nodespacing);
  reg.autoregister();

  reg.registered()->save_OIIO("registered.png");
  
}
