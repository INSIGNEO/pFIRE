#include "basewriter.hpp"

#include <string>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>

namespace bf = boost::filesystem;

const std::vector<std::string> BaseWriter::_components({"x", "y", "z"});

BaseWriter::BaseWriter(const std::string& filespec, const MPI_Comm &comm)
  : filename(BaseWriter::split_filespec(filespec).first),
    extra_path(BaseWriter::split_filespec(filespec).second), _comm(comm)
{
}

std::unique_ptr<BaseWriter::creator_map> BaseWriter::_creators =
    std::make_unique<BaseWriter::creator_map>();

std::unique_ptr<BaseWriter::extension_map> BaseWriter::_extension_handlers =
    std::make_unique<BaseWriter::extension_map>();

std::unique_ptr<stringset> BaseWriter::_truncated_files = std::make_unique<stringset>();

std::string BaseWriter::find_writer_for_filename(const std::string &filename)
{
  // Do it by lookup.
  string_pair pp = split_filespec(filename);
  std::string extension = bf::extension(split_filespec(filename).first);
  auto creatorname = BaseWriter::_extension_handlers->find(extension);
  if (creatorname != BaseWriter::_extension_handlers->end())
  {
    return creatorname->second;
  }
  std::ostringstream err;
  err << "No suitable writer available for extension \"" << extension << "\".";
  throw std::runtime_error(err.str());
}

BaseWriter_unique BaseWriter::get_writer_by_name(
    const std::string &name, const std::string &filename, MPI_Comm comm)
{
  // Do it by lookup.
  auto creator = BaseWriter::_creators->find(name);
  if (creator != BaseWriter::_creators->end())
  {
    return creator->second(filename, comm);
  }
  std::ostringstream err;
  err << "No registered writer named \"" << name << "\".";
  throw std::runtime_error(err.str());
}

BaseWriter_unique BaseWriter::get_writer_for_filename(const std::string &filename, MPI_Comm comm)
{
  try
  {
    return BaseWriter::get_writer_by_name(
        BaseWriter::find_writer_for_filename(filename), filename, comm);
  }
  catch (const std::runtime_error &errstr)
  {
    std::ostringstream err;
    err << "No suitable writer available for file \"" << filename << "\"."
        << "\n(" << errstr.what() << ")\n";
    throw std::runtime_error(err.str());
  }
}

BaseWriter::string_pair BaseWriter::split_filespec(std::string input)
{
  std::string::size_type delim_loc;
  delim_loc = input.find(':', 0);
  if (delim_loc != std::string::npos)
  {
    return std::make_pair(input.substr(0, delim_loc), input.substr(delim_loc + 1, input.size()));
  }

  return std::make_pair(input, std::string(""));
}

bool BaseWriter::check_truncated(const std::string& filename)
{
  if (!BaseWriter::_truncated_files)
  {
    _truncated_files = std::make_unique<stringset>();
  }
  return (_truncated_files->find(filename) != _truncated_files->end());
}

void BaseWriter::mark_truncated(std::string filename)
{
  if (!BaseWriter::_truncated_files)
  {
    _truncated_files = std::make_unique<stringset>();
  }
  _truncated_files->emplace(std::move(filename));
}
