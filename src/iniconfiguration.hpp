#ifndef INICONFIGURATION_HPP
#define INICONFIGURATION_HPP

#include <string>
#include <vector>

#include "baseconfiguration.hpp"
#include "types.hpp"

class IniConfig: public ConfigurationBase {
public:
  IniConfig(const int &argc, char const *const *argv);

  std::string usage();

  static bool valid_invocation(std::string &inv);

protected:
  void parse_arguments();

  void read_config_file(const std::string &config_path);
};
#endif // INICONFIGURATION_HPP
