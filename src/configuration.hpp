#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include <boost/filesystem.hpp>

namespace bf = boost::filesystem;

#include "types.hpp"

using config_map = std::map<std::string, std::string>;

class RegistrationConfig {
public:
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value, T>::type grab(const std::string key) const
  {
    return std::stoll(config.at(key));
  }

  template <typename T>
  typename std::enable_if<std::is_same<T, std::string>::value, T>::type grab(std::string key) const
  {
    return config.at(key);
  }

  static std::string get_invocation_name(const std::string& argzero);

  void validate_config();

protected:
  RegistrationConfig(const int& argc, char const* const* argv);

  config_map config;
  std::vector<std::string> arguments;
  std::string invocation_name;

  static const config_map default_config;
  static const std::vector<std::string> required_options;
  static const std::vector<std::string> arg_options;
  static const std::vector<std::string> bool_options;
};

#endif // CONFIGURATION_HPP
