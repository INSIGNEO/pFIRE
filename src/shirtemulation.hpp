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

#ifndef SHIRTEMULATION_HPP
#define SHIRTEMULATION_HPP

#include <string>
#include <vector>

#include "baseconfiguration.hpp"
#include "types.hpp"

class ShirtConfig: public ConfigurationBase {
public:
  ShirtConfig(const int& argc, char const* const* argv);

  std::string usage();

  static bool valid_invocation(std::string& inv);

protected:
  static const std::vector<std::string> shirt_synonyms;

  void parse_arguments();
};
#endif // SHIRTEMULATION_HPP
