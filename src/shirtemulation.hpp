#ifndef SHIRTEMULATION_HPP
#define SHIRTEMULATION_HPP

#include<string>
#include<vector>

#include "types.hpp"
#include "configuration.hpp"

class ShirtConfig: public RegistrationConfig
{
public:
  ShirtConfig(const int &argc, char const* const* argv);

  std::string usage();

  static bool valid_invocation(std::string &inv);

protected:
    
  static const std::vector<std::string> shirt_synonyms;

  void parse_arguments();

};
#endif //SHIRTEMULATION_HPP
