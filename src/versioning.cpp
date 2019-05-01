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

#include <petscsys.h>
#include <boost/version.hpp>

#ifdef USE_DCMTK
#include <dcmtk/config/osconfig.h>
#endif //USE_DCMTK

#ifdef USE_OIIO
#include <OpenImageIO/oiioversion.h>
#endif //USE_OIIO

#include "types.hpp"

std::string get_version_string_petsc()
{
#if PETSC_VERSION_MINOR > 7
    integer major, minor, subminor;
    PetscGetVersionNumber(&major, &minor, &subminor, nullptr);
        
    return std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(subminor);
#else
    return std::to_string(PETSC_VERSION_MAJOR) + "." + std::to_string(PETSC_VERSION_MINOR)
      + "." + std::to_string(PETSC_VERSION_SUBMINOR);
#endif

}

std::string get_version_string_boost()
{
  // See "Boost Informational Macros" section of boost::config docs
  integer major = BOOST_VERSION / 100000;
  integer minor = (BOOST_VERSION / 100) % 1000;
  integer subminor = BOOST_VERSION % 100;

  return std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(subminor);
}

#ifdef USE_DCMTK
std::string get_version_string_dcmtk()
{
  return std::string(PACKAGE_VERSION);
}
#endif //USE_DCMTK

#ifdef USE_OIIO
std::string get_version_string_oiio()
{
  return std::to_string(OIIO_VERSION_MAJOR) + "." + std::to_string(OIIO_VERSION_MINOR) + "." 
    + std::to_string(OIIO_VERSION_PATCH);
}
#endif //USE_OIIO
