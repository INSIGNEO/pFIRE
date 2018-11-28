#include "shirtemulation.hpp"

#include<cstdlib>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>

#include<boost/algorithm/string/case_conv.hpp>

namespace ba = boost::algorithm;

const std::vector<std::string> ShirtConfig::shirt_synonyms {"shirt", "pfire-shirt"};

std::string ShirtConfig::usage()
{
  std::ostringstream usagess;
  usagess << "Usage: "
          << invocation_name << " Register Fixed <filename> Moved <filename> Mask <filename> "
          << "NodeSpacing <nodespacing> [Map <filename>] [Registered <filename>]";

  return usagess.str();
}

ShirtConfig::ShirtConfig(const int &argc, char const* const* argv)
  : RegistrationConfig(argc, argv)
{
  // parent constructor provides arguments list so just need to parse
  parse_arguments();
}

void ShirtConfig::parse_arguments()
{
  if(arguments.size() == 0)
  {
    std::cout << usage() << std::endl;
    std::exit(0);
  }

  auto args_it = arguments.cbegin();
  // check mode is register or help
  std::string shirt_mode = ba::to_lower_copy(*args_it);
  if(shirt_mode == "help")
  {
    std::cout << usage() << std::endl;
    std::exit(0);
  }
  else if(shirt_mode != "register")
  {
    std::cout << "Unknown mode \"" << *args_it << "\".\n\n"
              << usage() << std::endl;
    std::exit(0);
  }
  std::advance(args_it, 1);
  while(args_it != arguments.cend())
  {
    //get argument and increment iterator
    std::string arg_lower = ba::to_lower_copy(*args_it);
    std::advance(args_it, 1); 
    
    //first check for "help"
    if(arg_lower == "help")
    {
      std::cout << usage() << std::endl;
      std::exit(0);
    }
    
    //now check option taking arguments
    if(std::find(arg_options.cbegin(), arg_options.cend(), arg_lower) != 
       arg_options.cend())
    {
      //try to get associated option
      if(args_it == arguments.cend())
      {
        //no argument, this is an error
        throw std::runtime_error("missing option for argument");
      }
      //grab associated option and increment iterator
      std::string optval = *args_it;
      std::advance(args_it, 1); 
      //finally insert arg-val pair into the argument map
      config[arg_lower] = optval;
    }
    else if(std::find(bool_options.cbegin(), bool_options.cend(), arg_lower) != 
            bool_options.cend())
    {
      //insert arg-val pair into argument map
      config[arg_lower] = "true";
    }
    else
      //if anything left over we have unhandled args
    {
      std::ostringstream errbuf;
      errbuf << "Unhandled arguments: ";
      std::copy(arguments.cbegin(), arguments.cend(), std::ostream_iterator<std::string>(errbuf, " "));
      throw std::runtime_error(errbuf.str());
    }
  }
}

bool ShirtConfig::valid_invocation(std::string &inv)
{
  if(std::find(shirt_synonyms.cbegin(), shirt_synonyms.cend(), inv) != shirt_synonyms.end())
  {
    return true;
  }
  return false;
}


