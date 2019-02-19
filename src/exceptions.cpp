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

#include <iostream>

#include "gitstate.hpp"

constexpr std::string_view abort_info_pre = R"FATAL(!!! FATAL ERROR !!!:

An unhandled exception has occured within the program.  This should not have happened and is 
probably a bug.

Help is available through the pFIRE issue tracker on github, but first please ensure:

1. You are using the most recent pFIRE release (check https://github.com/INSIGNEO/pFIRE/releases)
2. Your input files are not corrupt
3. If you are running on HPC please check there are no hardware issues

If you have checked the above and the bug still happens: https://github.com/INSIGNEO/pFIRE/issues

Your report should include:

- Input configuration file
- Details of the images being registered
- The exact pFIRE version: )FATAL";

constexpr std::string_view abort_info_post = R"FATAL(


)FATAL";


void abort_with_unhandled_error(){
  
  std::cout << abort_info_pre
            << git_version_string()
            << abort_info_post;

  std::abort();

}

