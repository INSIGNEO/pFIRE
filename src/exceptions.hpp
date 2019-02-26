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

#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <sstream>

#include "types.hpp"

void abort_with_unhandled_error();
void sigterm_handler(int signal);

class InvalidLoaderError: public std::runtime_error {
public:
  InvalidLoaderError(const std::string& filepath) : std::runtime_error(build_errstring(filepath))
  {
  }

protected:
  static std::string build_errstring(const std::string& path)
  {
    std::ostringstream errss;
    errss << "Failed to read data from " << path << ", wrong loader or file corrupt.";
    return errss.str();
  }
};

class FileNotFoundError: public std::runtime_error {
public:
  FileNotFoundError(const std::string& filepath) : std::runtime_error(build_errstring(filepath)) {}

protected:
  static std::string build_errstring(const std::string& path)
  {
    std::ostringstream errss;
    errss << "Failed to open " << path << ", wrong permissions or file does not exist.";
    return errss.str();
  }
};

class InternalError: public std::runtime_error {
public:
  InternalError(const std::string& what, std::string file = "unknown", integer line = 0)
    : std::runtime_error(build_errstring(what, file, line))
  {
  }

protected:
  static std::string build_errstring(
      const std::string& what, std::string file, integer line)
  {
    std::ostringstream errss;
    errss << "Internal error at " << file << ":" << line << " \"" << what << "\"";
    return errss.str();
  }
};

#endif // EXCEPTIONS_HPP
