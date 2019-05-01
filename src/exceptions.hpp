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
void print_abort_message();

class pFIREBaseError: public std::runtime_error {
public:
  pFIREBaseError(const std::string &msg, const std::string &file, size_t line_num)
    : std::runtime_error(msg), _where(build_wherestring(file, line_num))
  {
  }

  const char *where() { return _where.c_str(); }

protected:
  std::string _where;

  static std::string build_wherestring(const std::string &file, size_t line_num)
  {
    std::ostringstream errss;
    errss << file << ":" << line_num;
    return errss.str();
  }

  static std::string build_errstring(const std::string &what)
  {
    std::ostringstream errss;
    errss << what;
    return errss.str();
  }

};

class pFIREExpectedError: public pFIREBaseError {
public:
  pFIREExpectedError(const std::string &msg, const std::string& file, size_t line_num)
    : pFIREBaseError(msg, file, line_num)
  {}

};

class InvalidLoaderError: public pFIREExpectedError {
public:
  InvalidLoaderError(const std::string &path, std::string file = "unknown", size_t line_num = 0)
    : pFIREExpectedError(build_errstring(path), file, line_num)
  {
  }

protected:
  static std::string build_errstring(const std::string &path)
  {
    std::ostringstream errss;
    errss << "Failed to read data from " << path << ", wrong loader or file corrupt.";
    return errss.str();
  }
};

class FileNotFoundError: public pFIREExpectedError {
public:
  FileNotFoundError(const std::string &path, std::string file = "unknown", size_t line_num = 0)
    : pFIREExpectedError(build_errstring(path), file, line_num)
  {
  }

protected:
  static std::string build_errstring(const std::string &path)
  {
    std::ostringstream errss;
    errss << "Failed to open " << path << ", wrong permissions or file does not exist.";
    return errss.str();
  }
};

class InternalError: public pFIREBaseError {
public:
  InternalError(const std::string &what, std::string file = "unknown", size_t line = 0)
    : pFIREBaseError(build_errstring(what, file, line), file, line)
  {
  }

protected:
  std::string build_errstring(const std::string &what, const std::string &file, size_t line)
  {
    std::ostringstream errss;
    errss << "Internal error at " << file << ":" << line << " \"" << what << "\"";
    return errss.str();
  }
};

class BadConfigurationError: public pFIREExpectedError {
public:
  BadConfigurationError(const std::string &msg, std::string file = "unknown", size_t line = 0)
    : pFIREExpectedError(build_errstring(msg), file, line)
  {
  }

};

class InvalidWriterError: public pFIREExpectedError {
public:
  InvalidWriterError(const std::string &msg, std::string file = "unknown", size_t line = 0)
    : pFIREExpectedError(build_errstring(msg), file, line)
  {
  }
};


class WriterError: public pFIREExpectedError {
public:
  WriterError(const std::string &msg, std::string file = "unknown", size_t line = 0)
    : pFIREExpectedError(build_errstring(msg), file, line)
  {
  }
};

#endif // EXCEPTIONS_HPP
