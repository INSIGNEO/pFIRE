#include "basewriter.hpp"

#include <string>
#include <vector>

const std::vector<std::string> BaseWriter::_components({"x", "y", "z"});

BaseWriter::BaseWriter(std::string filename, const MPI_Comm& comm)
  : filename(std::move(filename)), _comm(comm)
{
}
