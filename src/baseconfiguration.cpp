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

#include "baseconfiguration.hpp"

#include <sstream>

#include "exceptions.hpp"
#include "infix_iterator.hpp"

const std::string ConfigurationBase::k_stem_token = "%name%";
const std::string ConfigurationBase::k_extension_token = "%ext%";
const std::string ConfigurationBase::k_outer_token = "%s%";
const std::string ConfigurationBase::k_inner_token = "%i%";

const config_map ConfigurationBase::arg_options = {
    {"registered", "registered.xdmf"}, {"map", "map.xdmf"}, {"lambda", "auto"}, {"mask", ""},
    {"registered_h5_path", "/registered"}, {"map_h5_path", "/map"}, {"lambda_mult", "1.0"},
    {"intermediate_template", "%name%-intermediate-%s%-%i%%ext%"},
    {"intermediate_map_template", "%name%-intermediate-map-%s%-%i%%ext%"},
    {"intermediate_directory", "intermediates"}, {"max_iterations", "100"}};

const config_map ConfigurationBase::bool_options = {
  {"verbose", "false"}, {"with_memory", "true"}, {"save_intermediate_frames", "false"}};

const std::vector<std::string> ConfigurationBase::required_options = {
    "fixed", "moved", "nodespacing"};

ConfigurationBase::ConfigurationBase(const int &argc, char const *const *argv)
  : arguments(argv + 1, argv + argc), invocation_name(get_invocation_name(argv[0]))
{
  config.insert(arg_options.cbegin(), arg_options.cend());
  config.insert(bool_options.cbegin(), bool_options.cend());
}

void ConfigurationBase::validate_config()
{
  std::list<std::string> missing;
  for (auto &req_it : ConfigurationBase::required_options)
  {
    if (config.find(req_it) == config.cend())
    {
      missing.push_back(req_it);
    }
  }
  if (missing.empty())
  {
    std::ostringstream errmsg;
    errmsg << "Missing required argument(s) \"";
    std::copy(missing.cbegin(), missing.cend(), infix_ostream_iterator<std::string>(errmsg, ", "));
    errmsg << "\"";
    throw BadConfigurationError(errmsg.str());
  }
}

std::string ConfigurationBase::get_invocation_name(const std::string &argzero)
{
  bf::path invpath(argzero);
  return (invpath.filename().string());
}
