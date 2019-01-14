#include "basewriter.hpp"

#include <boost/filesystem.hpp>

#include <string>
#include <vector>

namespace bf = boost::filesystem;

const std::vector<std::string> BaseWriter::_components({"x", "y", "z"});

BaseWriter::BaseWriter(std::string filename, const MPI_Comm &comm)
  : filename(std::move(filename)), _comm(comm)
{
}

std::unique_ptr<BaseWriter::creator_map> BaseWriter::_creators =
    std::make_unique<BaseWriter::creator_map>();

std::unique_ptr<BaseWriter::extension_map> BaseWriter::_extension_handlers =
    std::make_unique<BaseWriter::extension_map>();

std::string BaseWriter::find_writer_for_filename(const std::string &filename)
{
  // Do it by lookup.
  std::string extension = bf::extension(filename);
  auto creatorname = BaseWriter::_extension_handlers->find(extension);
  if (creatorname != BaseWriter::_extension_handlers->end())
  {
    return creatorname->second;
  }
  std::ostringstream err;
  err << "No suitable writer available for extension \"" << extension << "\".";
  throw std::runtime_error(err.str());
}

BaseWriter_unique BaseWriter::get_writer_by_name(const std::string &name,
                                                 const std::string &filename, MPI_Comm comm)
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
  // Do it by lookup.
  std::string extension = bf::extension(filename);
  try
  {
    return BaseWriter::get_writer_by_name(BaseWriter::find_writer_for_filename(filename),
                                          filename, comm);
  }
  catch (const std::runtime_error&) 
  {
  std::ostringstream err;
  err << "No suitable writer available for extension \"" << extension << "\".";
  throw std::runtime_error(err.str());
  }
}
