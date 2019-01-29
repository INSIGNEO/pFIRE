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

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <boost/filesystem.hpp>

namespace bf = boost::filesystem;

#include "types.hpp"

using config_map = std::map<std::string, std::string>;

class ConfigurationBase {
public:

  template <typename T>
  typename std::enable_if<std::is_same<T, bool>::value, T>::type grab(const std::string key) const
  {
    T val;
    std::istringstream(config.at(key)) >> std::boolalpha >> val;
    return val;
  }

  template <typename T>
  typename std::enable_if<std::is_integral<T>::value && !std::is_same<T, bool>::value, T>::type
  grab(const std::string key) const
  {
    return static_cast<T>(std::stoll(config.at(key)));
  }

  template <typename T>
  typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, bool>::value, T>::type
  grab(const std::string key) const
  {
    return static_cast<T>(std::stold(config.at(key)));
  }

  template <typename T>
  typename std::enable_if<std::is_same<T, std::string>::value, T>::type grab(std::string key) const
  {
    return config.at(key);
  }

  static std::string get_invocation_name(const std::string& argzero);

  void validate_config();

protected:
  ConfigurationBase(const int& argc, char const* const* argv);

  config_map config;
  std::vector<std::string> arguments;
  std::string invocation_name;

  static const config_map default_config;
  static const std::vector<std::string> required_options;
  static const std::vector<std::string> arg_options;
  static const std::vector<std::string> bool_options;
};

#endif // CONFIGURATION_HPP
