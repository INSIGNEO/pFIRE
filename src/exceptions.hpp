#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <sstream>

class InvalidLoaderError : public std::runtime_error {
public:
  InvalidLoaderError(const std::string& filepath)
  : std::runtime_error(file_to_string(filepath)) {}

protected:
  static std::string file_to_string(const std::string& path)
  {
    std::ostringstream errss;
    errss << "Failed to read data from " << path << ", wrong loader or file corrupt.";
    return errss.str();
  }
};

class FileNotFoundError : public std::runtime_error {
public:
  FileNotFoundError(const std::string& filepath)
  : std::runtime_error(file_to_string(filepath)) {}

protected:
  static std::string file_to_string(const std::string& path)
  {
    std::ostringstream errss;
    errss << "Failed to open " << path << ", wrong permissions or file does not exist.";
    return errss.str();
  }
};

#endif // EXCEPTIONS_HPP
