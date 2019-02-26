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

#include "setup.hpp"

#include <csignal>

#include <petscsys.h>

#include <boost/filesystem.hpp>

#include "types.hpp"
#include "banner.hpp"
#include "gitstate.hpp"
#include "baseloader.hpp"
#include "shirtloader.hpp"

#include "exceptions.hpp"

#ifdef USE_OIIO
#include "oiioloader.hpp"
#include "oiiowriter.hpp"
#endif // USE_OIIO

#ifdef USE_DCMTK
#include "dcmloader.hpp"
#endif // USE_DCMTK

#include "basewriter.hpp"
#include "hdfwriter.hpp"
#include "xdmfwriter.hpp"

namespace bf = boost::filesystem;

void register_plugins()
{
#ifdef USE_DCMTK
  BaseLoader::register_loader(DCMLoader::loader_name, DCMLoader::Create_Loader);
#endif // USE_DCMTK

#ifdef USE_OIIO
  BaseLoader::register_loader(OIIOLoader::loader_name, OIIOLoader::Create_Loader);
  BaseWriter::register_writer<OIIOWriter>();
#endif // USE_OIIO

  BaseLoader::register_loader(ShIRTLoader::loader_name, ShIRTLoader::Create_Loader);

  BaseWriter::register_writer<HDFWriter>();
  BaseWriter::register_writer<XDMFWriter>();
}

void pfire_setup(const std::vector<std::string>& petsc_args)
{
  // Setup terminate handler
  std::set_terminate(abort_with_unhandled_error);
  std::signal(SIGTERM, sigterm_handler);

  std::vector<char*> cstrings;
  cstrings.resize(petsc_args.size());
  std::transform(
      petsc_args.cbegin(), petsc_args.cend(), cstrings.begin(),
      [](auto str_it) -> char* { return const_cast<char*>(str_it.c_str()); });
  int n_cstrings = cstrings.size();
  char** cstr_ptr = cstrings.data();

  PetscErrorCode perr = PetscInitialize(&n_cstrings, &cstr_ptr, nullptr, nullptr);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(rank == 0)
  {
    print_welcome_message();
  }

  check_comm_size_and_warn_odd();

  register_plugins();
}

void check_comm_size_and_warn_odd()
{
  // Check for even number of processors, warn on oddness
  int comm_size, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(rank == 0)
  {
    std::cout << "Running on " << comm_size << " processes" << std::endl;

    if (comm_size > 2 && comm_size % 2 == 1)
    {
      std::cout << "!!!! WARNING !!!!\n"
                << "Using an odd number of processors is not recommended.\n"
                << "This makes efficient subdivision of the problem much harder and will likely "
                << "lead to reduced performance.\n"
                << "Extreme cases may make viable partitioning impossible and cause job failure. "
                << "You have been warned!\n\n"
                << std::flush;
    }
  }
}

void pfire_teardown()
{
  PetscFinalize();
}


void print_welcome_message()
{

  std::ostringstream welcomess;
  welcomess << kbanner_text_upper;

  welcomess << git_version_string();

  welcomess << kbanner_text_lower;

  std::cout << welcomess.str();
}
