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

#include "baseconfiguration.hpp"
#include "iniconfiguration.hpp"
#include "setup.hpp"
#include "shirtemulation.hpp"

#include "elastic.hpp"
#include "image.hpp"
#include "infix_iterator.hpp"
#include "laplacian.hpp"
#include "map.hpp"
#include "types.hpp"
#include "math_utils.hpp"

#include "xdmfwriter.hpp"

void mainflow(std::shared_ptr<ConfigurationBase> config);

void usage()
{
  PetscPrintf(PETSC_COMM_WORLD, "Usage: pfire fixed moved nodespacing\n");
}

int main(int argc, char **argv)
{
  pfire_setup(std::vector<std::string>());

  std::string invocation_name = ConfigurationBase::get_invocation_name(argv[0]);

  std::shared_ptr<ConfigurationBase> configobj(nullptr);
  if (ShirtConfig::valid_invocation(invocation_name))
  {
    configobj = std::make_shared<ShirtConfig>(argc, argv);
  }
  else
  {
    configobj = std::make_shared<IniConfig>(argc, argv);
  }

  configobj->validate_config();

  auto tstart = std::chrono::high_resolution_clock::now();
  mainflow(configobj);

  auto tend = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff = tend - tstart;

  PetscPrintf(PETSC_COMM_WORLD, "Elapsed time: %g s\n", diff.count());

  pfire_teardown();

  return 0;
}

void mainflow(std::shared_ptr<ConfigurationBase> config)
{
  std::unique_ptr<Image> fixed;
  try
  {
    fixed = Image::load_file(config->grab<std::string>("fixed"));
  }
  catch (std::exception &e)
  {
    std::cerr << "Error: Failed to load fixed image: " << e.what() << std::endl;
    return;
  }

  std::ostringstream immsg;
  immsg << "Loaded fixed image of shape ";
  std::copy_n(
      fixed->shape().cbegin(), fixed->ndim(), infix_ostream_iterator<integer>(immsg, " x "));
  immsg << ".\n";
  PetscPrintf(PETSC_COMM_WORLD, immsg.str().c_str());

  std::unique_ptr<Image> moved;
  try
  {
    moved = Image::load_file(config->grab<std::string>("moved"), fixed.get());
  }
  catch (std::exception &e)
  {
    std::cerr << "Error: Failed to load moved image: " << e.what() << std::endl;
    return;
  }

  floatvector nodespacing(fixed->ndim(), config->grab<integer>("nodespacing"));

  // explain_memory(fixed->shape(), Map::calculate_map_shape(fixed->shape(), nodespacing));

  fixed->normalize();
  moved->normalize();

  std::unique_ptr<Mask> mask;
  if(config->grab<std::string>("mask").empty()){
    mask = Mask::full_image(*fixed);
  }
  else
  {
    try
    {
      mask = Mask::load_file(config->grab<std::string>("mask"), fixed.get());
    }
    catch (std::exception &e)
    {
      std::cerr << "Error: Failed to load mask: " << e.what() << std::endl;
      return;
    }
  }

  Elastic reg(*fixed, *moved, *mask, nodespacing, *config);
  reg.autoregister();

  std::string outfile = config->grab<std::string>("registered");
  BaseWriter_unique wtr = BaseWriter::get_writer_for_filename(outfile, fixed->comm());
  wtr->write_image(*reg.registered());

  outfile = config->grab<std::string>("map");
  wtr = BaseWriter::get_writer_for_filename(outfile, fixed->comm());
  wtr->write_map(*reg.m_p_map);
}
