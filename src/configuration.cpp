#include "configuration.hpp"

#include<sstream>

#include "infix_iterator.hpp"

const config_map RegistrationConfig::default_config = {
  {"verbose", "false"},
  {"registered", "registered"},
  {"map", "map"}};

const std::vector<std::string> RegistrationConfig::required_options = {
  "fixed",
  "moved",
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
  std::list<std::string> missing;
  for(auto &req_it: RegistrationConfig::required_options)
  {
    if(config.find(req_it) == config.cend())
    {
      missing.push_back(req_it);
    }
  }
  if(missing.size() > 0)
  {
    std::ostringstream errmsg;
    errmsg << "Missing required argument(s) \"";
    std::copy(missing.cbegin(), missing.cend(), infix_ostream_iterator<std::string>(errmsg, ", "));
    errmsg << "\"";
    throw std::runtime_error(errmsg.str());
  }
}

std::string RegistrationConfig::get_invocation_name(const std::string &argzero)
{
  bf::path invpath(argzero);
  return(invpath.filename().string());
}
