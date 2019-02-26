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

#ifndef SETUP_HPP
#define SETUP_HPP

#include <string>
#include <vector>

std::string get_invocation_name(const std::string &argzero);

void pfire_setup(const std::vector<std::string> &petsc_args);

void check_comm_size_and_warn_odd();
void print_welcome_message();

void pfire_teardown();

#endif // SETUP_HPP
