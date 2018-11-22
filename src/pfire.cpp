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

#include "baseloader.hpp"

#ifdef USE_DCMTK
#include "dcmloader.hpp"
#endif //USE_DCMTK

#ifdef USE_OIIO
#include "oiioloader.hpp"
#endif //USE_OIIO

void mainflow(std::string, std::string, floating);

void register_plugins()
{
#ifdef USE_DCMTK
  BaseLoader::register_loader(DCMLoader::loader_name, DCMLoader::Create_Loader);
#endif //USE_DCMTK

#ifdef USE_OIIO
  BaseLoader::register_loader(OIIOLoader::loader_name, OIIOLoader::Create_Loader);
#endif //USE_OIIO
}

void usage()
{
  PetscPrintf(PETSC_COMM_WORLD, "Usage: pfire fixed moved nodespacing\n");
}

int main(int argc, char **argv){
  PetscErrorCode perr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRXX(perr);

  register_plugins();

  if(argc != 4)
  {
    usage();
    return 0;
  }

  floating nodespacing = std::stod(argv[3]);

  mainflow(argv[1], argv[2], nodespacing);

  PetscFinalize();

  return 0;
}

void mainflow(std::string fixedpath, std::string movedpath, floating ns){


  auto foo = BaseLoader::loaders();


  std::unique_ptr<Image> fixed;
  try
  {
    fixed = Image::load_file(fixedpath);
  }
  catch(std::exception &e)
  {
    std::cerr << "Error: Failed to load fixed image: " << e.what() << std::endl;
    return;
  }
  PetscPrintf(PETSC_COMM_WORLD, "Loaded fixed image of shape %i x %i x %i.\n", fixed->shape()[0],
              fixed->shape()[1], fixed->shape()[2]);

  std::unique_ptr<Image> moved;
  try
  {
    moved = Image::load_file(movedpath);
  }
  catch(std::exception &e)
  {
    std::cerr << "Error: Failed to load moved image: " << e.what() << std::endl;
    return;
  }

  fixed->normalize();
  moved->normalize();

  floatvector nodespacing(fixed->ndim(), ns);

  Elastic reg(*fixed, *moved, nodespacing);
  reg.autoregister();

}
