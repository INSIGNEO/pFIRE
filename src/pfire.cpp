#include <chrono>

#include "setup.hpp"
#include "configuration.hpp"
#include "shirtemulation.hpp"
#include "iniconfiguration.hpp"

#include "types.hpp"
#include "laplacian.hpp"
#include "image.hpp"
#include "map.hpp"
#include "elastic.hpp"
#include "utils.hpp"

#include "hdfwriter.hpp"

void mainflow(std::string, std::string, floating);

void usage()
{
  PetscPrintf(PETSC_COMM_WORLD, "Usage: pfire fixed moved nodespacing\n");
}

int main(int argc, char **argv){

  pfire_setup(std::vector<std::string>());

  std::string invocation_name = RegistrationConfig::get_invocation_name(argv[0]);

  std::unique_ptr<RegistrationConfig> configobj(nullptr);

  if(ShirtConfig::valid_invocation(invocation_name))
  {
    configobj = std::make_unique<ShirtConfig>(argc, argv);
  }
  else
  {
    configobj = std::make_unique<IniConfig>(argc, argv);
  }

  configobj->validate_config();

  auto tstart = std::chrono::high_resolution_clock::now();
  mainflow(configobj->grab<std::string>("fixed"),
           configobj->grab<std::string>("moved"),
           configobj->grab<integer>("nodespacing"));

  auto tend = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff = tend - tstart;

  PetscPrintf(PETSC_COMM_WORLD, "Elapsed time: %g s\n", diff.count());

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

  floatvector nodespacing(fixed->ndim(), ns);

  explain_memory(fixed->shape(), Map::calculate_map_shape(fixed->shape(), nodespacing));

  fixed->normalize();
  moved->normalize();


  Elastic reg(*fixed, *moved, nodespacing);
  reg.autoregister();

  HDFWriter wtr("data.h5", fixed->comm());

  std::string reggroup("registered");
  std::string mapgroup("map");
  wtr.write_image(*reg.registered(), reggroup);
  wtr.write_map(*reg.m_p_map, mapgroup);
  
}
