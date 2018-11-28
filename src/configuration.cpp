#include "configuration.hpp"

#include<sstream>

const config_map RegistrationConfig::default_config = {
  {"verbose", "false"},
  {"registered", "registered"},
  {"map", "map"}};

const std::vector<std::string> RegistrationConfig::required_options = {
  "fixed",
  "moved",
  "mask",
  "nodespacing"};

const std::vector<std::string> RegistrationConfig::arg_options = {
  "fixed",
  "moved",
  "mask", 
  "nodespacing"};

const std::vector<std::string> RegistrationConfig::bool_options = {
  "verbose"};

RegistrationConfig::RegistrationConfig(const int &argc, char const* const* argv)
  : config(default_config), arguments(argv+1, argv+argc)
{
  invocation_name = get_invocation_name(argv[0]);
}

void RegistrationConfig::validate_config()
{
  for(auto &req_it: RegistrationConfig::required_options)
  {
    if(config.find(req_it) == config.cend())
    {
      std::ostringstream errmsg;
      errmsg << "Missing required argument " << req_it;
      throw std::runtime_error(errmsg.str());
    }
  }
}

std::string RegistrationConfig::get_invocation_name(const std::string &argzero)
{
  bf::path invpath(argzero);
  return(invpath.filename().string());
}
