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

#include "exceptions.hpp"

#include <chrono>
#include <iostream>
#include <thread>

#include <mpi.h>

#include "gitstate.hpp"
#include "versioning.hpp"

constexpr std::string_view abort_info_pre = R"FATAL(!!! FATAL ERROR !!!:

An unhandled exception has occured within the program.  This should not have happened and is 
probably a bug.

Help is available through the pFIRE issue tracker on github, but first please ensure:

1. You are using the most recent pFIRE release (check https://github.com/INSIGNEO/pFIRE/releases)
2. Your input files are not corrupt
3. If you are running on HPC please check there are no hardware issues

If you have checked the above and the bug still occurs please report it via 
https://github.com/INSIGNEO/pFIRE/issues

Your report should include:

- Input configuration file
- Details of the images being registered (format, dimensions)
- Basic hardware details, operating system and version
- The exact pFIRE version:

    - )FATAL";

constexpr std::string_view abort_info_libheader = R"FATAL(

- The library versions in use:

)FATAL";

constexpr std::string_view abort_info_petsc = "    - Petsc Version: ";
constexpr std::string_view abort_info_boost = "    - Boost Version: ";

#ifdef USE_OIIO
constexpr std::string_view abort_info_oiio = "    - OpenImageIO Version: ";
#endif

#ifdef USE_DCMTK
constexpr std::string_view abort_info_dcmtk = "    - DCMTK Version: ";
#endif //USE_DCMTK

constexpr std::string_view abort_info_post = R"FATAL(
)FATAL";


void abort_with_unhandled_error(){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::exception_ptr eptr = std::current_exception();

  // Try to separate threads that arrive here together
  // to avoid garbling stdout too much
  const int waittime(200);
  std::this_thread::sleep_for(std::chrono::milliseconds(waittime*rank));

  try
  {
    if (eptr) {
      std::rethrow_exception(eptr);
    }
  }
  catch(const std::exception& e)
  {
    // Use std::cout here because sensible MPI behaviour has gone out the window
    std::cout << "Rank " << rank << " encountered an unhandled exception: \n\"\"\"\n" 
              << e.what() << "\n\"\"\"\n";
  }

  print_abort_message();

  MPI_Abort(MPI_COMM_WORLD, -1); 

}

void sigterm_handler(int signal){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Use std::cout here because sensible MPI behaviour has gone out the window
  std::cout << "Rank " << rank << " received signal " << signal;

  const int waittime(200);
  std::this_thread::sleep_for(std::chrono::milliseconds(waittime*rank));
  print_abort_message();
}

void print_abort_message(){

  // Use std::cout here because sensible MPI behaviour has gone out the window
  std::cout << abort_info_pre
            << git_version_string()
            << abort_info_libheader
            << abort_info_petsc << get_version_string_petsc() << "\n"
            << abort_info_boost << get_version_string_boost() << "\n"
#ifdef USE_OIIO
            << abort_info_oiio << get_version_string_oiio() << "\n"
#endif //USE_OIIO
#ifdef USE_DCMTK
            << abort_info_dcmtk << get_version_string_dcmtk() << "\n"
#endif //USE_DCMTK
            << abort_info_post;
}
