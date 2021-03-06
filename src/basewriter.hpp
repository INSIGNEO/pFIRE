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

#ifndef BASEWRITER_HPP
#define BASEWRITER_HPP

#include <functional>
#include <map>
#include <string>
#include <type_traits>

#include <hdf5.h>
#include <mpi.h>

#include "types.hpp"
#include "exceptions.hpp"

class BaseWriter {
public:
  using writer_creator =
      std::function<BaseWriter_unique(const std::string& path, MPI_Comm comm)>;

  using creator_map = std::map<std::string, writer_creator>;
  using extension_map = std::map<std::string, std::string>;

  using string_pair = std::pair<std::string, std::string>;

  BaseWriter(const std::string& filespec, const MPI_Comm& comm = PETSC_COMM_WORLD);
  virtual ~BaseWriter() = default;

  virtual std::string write_image(const Image& image) = 0;
  virtual std::string write_map(const MapBase& map) = 0;

  template <typename WriterClass>
  static BaseWriter_unique create_writer(const std::string &path, MPI_Comm comm)
  {
    return BaseWriter_unique(new WriterClass(path, comm));
  }

  const std::string& filename() const { return _filename; }
  const std::string& extra_path() const { return _extra_path; }
  const MPI_Comm& comm() const { return _comm;}

  static const std::vector<std::string>& components() { return _components;}

  template <typename WriterClass>
  static bool
  //typename std::enable_if<std::is_convertible<WriterClass*, BaseWriter*>::value, bool>::type 
  register_writer();

  static std::string find_writer_for_filename(const std::string& _filename);

  static BaseWriter_unique get_writer_by_name(const std::string& name,
                                              const std::string& _filename,
                                              MPI_Comm comm = PETSC_COMM_WORLD);

  static BaseWriter_unique get_writer_for_filename(const std::string& _filename,
                                                    MPI_Comm comm = PETSC_COMM_WORLD);

  static string_pair split_filespec(std::string input);

  static bool check_truncated(const std::string& _filename);
  static void mark_truncated(std::string _filename);

  static const std::string writer_name;
  static const std::vector<std::string> extensions;

private:
  std::string _filename;
  std::string _extra_path;
  MPI_Comm _comm;

  static const std::vector<std::string> _components;

  static std::unique_ptr<creator_map> _creators;
  static std::unique_ptr<extension_map> _extension_handlers;
  static std::unique_ptr<stringset> _truncated_files;
};


template <typename WriterClass>
bool
//typename std::enable_if<std::is_convertible<WriterClass*, BaseWriter*>::value, bool>::type
BaseWriter::register_writer()
{
  if (!BaseWriter::_creators)
  {
    _creators = std::make_unique<BaseWriter::creator_map>();
  }
  if (!BaseWriter::_extension_handlers)
  {
    _extension_handlers = std::make_unique<BaseWriter::extension_map>();
  }
  if (auto it = _creators->find(WriterClass::writer_name); it == _creators->cend())
  {
    (*_creators)[WriterClass::writer_name] = BaseWriter::create_writer<WriterClass>;
    for (const auto &extit : WriterClass::extensions){
      BaseWriter::_extension_handlers->insert(extension_map::value_type(extit, WriterClass::writer_name));
    }
    return true;
  }
  std::ostringstream err;
  err << "Writer named " << WriterClass::writer_name << " already registered.";
  throw InternalError(err.str(), __FILE__, __LINE__);
}

#endif // BASEWRITER_HPP
