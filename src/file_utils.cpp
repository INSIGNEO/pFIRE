#include "file_utils.hpp"

#include <boost/filesystem.hpp>

#include "exceptions.hpp"

namespace bf = boost::filesystem;

void throw_if_nonexistent(const std::string& filepath)
{
  if (!bf::exists(filepath))
  {
    throw FileNotFoundError(filepath);
  }
}


