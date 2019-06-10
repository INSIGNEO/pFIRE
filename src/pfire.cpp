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

#include <chrono>
#include <cstdlib>

#include "types.hpp"
#include "baseconfiguration.hpp"
#include "iniconfiguration.hpp"
#include "setup.hpp"
#include "shirtemulation.hpp"

#include "mapmanager.hpp"
#include "exceptions.hpp"
#include "image.hpp"
#include "mapbase.hpp"
#include "mask.hpp"
#include "elasticregistration.hpp"
#include "basewriter.hpp"

void mainflow(const std::shared_ptr<ConfigurationBase>& config);

int main(int argc, char **argv)
{
  pfire_setup(std::vector<std::string>());

  std::string invocation_name = ConfigurationBase::get_invocation_name(argv[0]);

  std::shared_ptr<ConfigurationBase> configobj(nullptr);
  try{
    if (ShirtConfig::valid_invocation(invocation_name))
    {
      configobj = std::make_shared<ShirtConfig>(argc, argv);
    }
    else
    {
      configobj = std::make_shared<IniConfig>(argc, argv);
    }
  }
  catch (const BadConfigurationError &err)
  {
    PetscPrintf(MPI_COMM_WORLD, "Failed to parse configuration: %s", err.what());
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  try
  {
    configobj->validate_config();
  }
  catch (const BadConfigurationError &e)
  {
    std::cerr << "Error: failed to parse configuration: " << e.what() << std::endl;
    return -1;
  }

  auto tstart = std::chrono::high_resolution_clock::now();
  mainflow(configobj);

  auto tend = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff = tend - tstart;

  PetscPrintf(MPI_COMM_WORLD, "Elapsed time: %g s\n", diff.count());

  if(configobj->grab<bool>("memory_report"))
  {
    vmem_report();
  }

  pfire_teardown();

  return 0;
}

void mainflow(const std::shared_ptr<ConfigurationBase>& config)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::shared_ptr<Image> fixed;
  try
  {
    fixed = ImageBase::load_image<Image>(config->grab<std::string>("fixed"));
  }
  catch (const pFIREExpectedError &e)
  {
    std::cerr << "Error: Failed to load fixed image: " << e.what() << std::endl;
    return;
  }

  std::shared_ptr<Image> moved;
  try
  {
    moved = ImageBase::load_image<Image>(config->grab<std::string>("moved"), fixed.get());
  }
  catch (const pFIREExpectedError &e)
  {
    std::cerr << "Error: Failed to load moved image: " << e.what() << std::endl;
    return;
  }

  intcoord mapshape = {config->grab<integer>("nodespacing"),
                       config->grab<integer>("nodespacing"),
                       config->grab<integer>("nodespacing")};

  auto mask = Mask::create_filled_mask(*fixed);

  auto reg = ElasticRegistration(fixed, moved, mask, mapshape, *config);

  auto minmax = fixed->minmax();
  PetscPrintf(MPI_COMM_WORLD, "Fixed minmax: %f, %f\n", minmax.first, minmax.second); 

  reg.autoregister();

  try
  {
    std::string outfile = config->grab<std::string>("registered");
    std::string h5group = config->grab<std::string>("registered_h5_path");
    std::string output_path = outfile + ":" + h5group;
    BaseWriter_unique wtr = BaseWriter::get_writer_for_filename(output_path, fixed->comm());
    output_path = wtr->write_image(*reg.registered());
    PetscPrintf(PETSC_COMM_WORLD, "Saved registered image to %s\n", output_path.c_str());

    outfile = config->grab<std::string>("map");
    h5group = config->grab<std::string>("map_h5_path");
    output_path = outfile + ":" + h5group;
    wtr = BaseWriter::get_writer_for_filename(output_path, fixed->comm());
    output_path = wtr->write_map(reg.result_map());
    PetscPrintf(PETSC_COMM_WORLD, "Saved map to %s\n", output_path.c_str());

  }
  catch (const pFIREExpectedError &e)
  {
    std::cerr << "Error: Failed to save results: " << e.what() << std::endl;
    return;
  }

}
