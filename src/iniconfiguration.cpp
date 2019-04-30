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

#include "iniconfiguration.hpp"

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace ba = boost::algorithm;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

std::string IniConfig::usage()
{
  std::ostringstream usagess;
  usagess << "pFIRE (Parallel Framework for Image Registration)\n\n"
          << "Usage: " << invocation_name << " <config_file> [-h,--help]\n\n"
          << "Options:\n"
          << "  config_file\t\tconfiguration file to read options from";

  return usagess.str();
}

IniConfig::IniConfig(int const &argc, char const *const *argv) : ConfigurationBase(argc, argv)
{
  // parent constructor provides arguments list so just need to parse
  parse_arguments();
}

void IniConfig::parse_arguments()
{
  po::options_description cmdline_visible;
  cmdline_visible.add_options()("help,h", "print this message");

  po::options_description cmdline_hidden("Hidden positional options");
  cmdline_hidden.add_options()("config-file", po::value<std::string>(), "configuration ini file");

  po::positional_options_description positional;
  positional.add("config-file", 1);

  po::options_description cmdline;
  cmdline.add(cmdline_visible).add(cmdline_hidden);

  po::variables_map vm;
  try
  {
    po::store(
        po::command_line_parser(arguments).options(cmdline).positional(positional).run(), vm);
  }
  catch (const po::error &err)
  {
    std::cerr << "Error: " << err.what() << "\n\n"
              << usage() << std::endl
              << cmdline_visible << std::endl;
    std::exit(0);
  }

  if (vm.count("help") > 0 || vm.count("config-file") != 1 )
  {
    // TODO improve to iterate over options
    std::cout << usage() << std::endl << cmdline_visible << std::endl;
    std::exit(0);
  }

  try
  {
    po::notify(vm);
  }
  catch (const po::error &err)
  {
    std::cout << err.what() << std::endl;
  }

  read_config_file(vm["config-file"].as<std::string>());
}

void IniConfig::read_config_file(const std::string &config_path)
{
  pt::ptree config_data;
  try
  {
    pt::read_ini(config_path, config_data);
  }
  catch (const pt::ini_parser_error &err)
  {
    std::cerr << "Failed to read config file " << err.what() << std::endl;
    std::exit(0);
  }

  std::list<std::string> unknowns;
  for (const auto &it : config_data)
  {
    std::string key = ba::to_lower_copy(it.first);
    if (arg_options.find(key) != arg_options.cend())
    {
      config[key] = it.second.data();
    }
    else if (bool_options.find(key) != bool_options.cend())
    {
      config[key] = boost::to_lower_copy(it.second.data());
    }
    else
    {
      unknowns.push_back(key);
    }
  }

  if (!unknowns.empty())
  {
    std::ostringstream unkss;
    std::copy(unknowns.begin(), unknowns.end(), std::ostream_iterator<std::string>(unkss, ", "));
    std::cout << "Warning: unknown options \"" << unkss.str() << "\"\n";
  }
}

// Ignore unused inv because needed in other derived classes
bool IniConfig::valid_invocation(std::string &inv __attribute__((unused)))
{
  return true;
}
