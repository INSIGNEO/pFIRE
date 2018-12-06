#include "setup.hpp"

#include <petscsys.h>

#include <boost/filesystem.hpp>

#include "baseloader.hpp"
#include "shirtloader.hpp"
#include "types.hpp"

#ifdef USE_OIIO
#include "oiioloader.hpp"
#endif // USE_OIIO

#ifdef USE_DCMTK
#include "dcmloader.hpp"
#endif // USE_DCMTK

namespace bf = boost::filesystem;

void register_plugins()
{
#ifdef USE_DCMTK
  BaseLoader::register_loader(DCMLoader::loader_name, DCMLoader::Create_Loader);
#endif // USE_DCMTK

#ifdef USE_OIIO
  BaseLoader::register_loader(OIIOLoader::loader_name, OIIOLoader::Create_Loader);
#endif // USE_OIIO

  BaseLoader::register_loader(ShIRTLoader::loader_name, ShIRTLoader::Create_Loader);
}

void pfire_setup(const std::vector<std::string>& petsc_args)
{
  std::vector<char*> cstrings;
  cstrings.resize(petsc_args.size());
  std::transform(
      petsc_args.cbegin(), petsc_args.cend(), cstrings.begin(),
      [](auto str_it) -> char* { return const_cast<char*>(str_it.c_str()); });
  int n_cstrings = cstrings.size();
  char** cstr_ptr = cstrings.data();

  PetscErrorCode perr = PetscInitialize(&n_cstrings, &cstr_ptr, nullptr, nullptr);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  check_and_warn_odd_comm();

  register_plugins();
}

void check_and_warn_odd_comm()
{
  // Check for even number of processors, warn on oddness
  int comm_size, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (comm_size > 2 && comm_size % 2 == 1 && rank == 0)
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

void pfire_teardown()
{
  PetscFinalize();
}
