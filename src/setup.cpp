#include "setup.hpp"

#include<petscsys.h>

#include<boost/filesystem.hpp>

#include "types.hpp"
#include "baseloader.hpp"
#include "shirtloader.hpp"

#ifdef USE_OIIO
#include "oiioloader.hpp"
#endif //USE_OIIO

#ifdef USE_DCMTK
#include "dcmloader.hpp"
#endif //USE_DCMTK

namespace bf = boost::filesystem;

void register_plugins()
{
#ifdef USE_DCMTK
  BaseLoader::register_loader(DCMLoader::loader_name, DCMLoader::Create_Loader);
#endif //USE_DCMTK

#ifdef USE_OIIO
  BaseLoader::register_loader(OIIOLoader::loader_name, OIIOLoader::Create_Loader);
#endif //USE_OIIO

  BaseLoader::register_loader(ShIRTLoader::loader_name, ShIRTLoader::Create_Loader);
}

void pfire_setup(const std::vector<std::string> &petsc_args)
{
  std::vector<char* > cstrings;
  cstrings.resize(petsc_args.size());
  std::transform(petsc_args.cbegin(), petsc_args.cend(), cstrings.begin(),
                 [](auto str_it) -> char*{return const_cast<char*>(str_it.c_str());});
  int n_cstrings = cstrings.size();
  char ** cstr_ptr = cstrings.data();

  PetscErrorCode perr = PetscInitialize(&n_cstrings, &cstr_ptr, nullptr, nullptr);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  
  register_plugins();
}

void pfire_teardown(){
  PetscFinalize(); 
}
