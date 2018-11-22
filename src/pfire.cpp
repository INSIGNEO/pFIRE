#include "setup.hpp"

#include "types.hpp"
#include "laplacian.hpp"
#include "image.hpp"
#include "map.hpp"
#include "elastic.hpp"

//#include "hdfwriter.hpp"

void mainflow(std::string, std::string, floating);

void usage()
{
  PetscPrintf(PETSC_COMM_WORLD, "Usage: pfire fixed moved nodespacing\n");
}

int main(int argc, char **argv){

  std::cout << get_invocation_name(argv[0]) << std::endl;

  pfire_setup(std::vector<std::string>());

  if(argc < 4)
  {
    usage();
    return 0;
  }

  floating nodespacing = std::stod(argv[3]);

  mainflow(argv[1], argv[2], nodespacing);

  pfire_teardown();

  return 0;
}

void mainflow(std::string fixedpath, std::string movedpath, floating ns){

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
    moved = Image::load_file(movedpath, fixed.get());
  }
  catch(std::exception &e)
  {
    std::cerr << "Error: Failed to load moved image: " << e.what() << std::endl;
    return;
  }

  fixed->save_OIIO("fixed_writeback.png");
  moved->save_OIIO("moved_writeback.png");

  fixed->normalize();
  moved->normalize();

  floatvector nodespacing(fixed->ndim(), ns);

  Elastic reg(*fixed, *moved, nodespacing);
  reg.autoregister();

  reg.registered()->save_OIIO("registered.png");
  
}
