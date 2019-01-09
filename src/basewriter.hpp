#ifndef BASEWRITER_HPP
#define BASEWRITER_HPP

#include <string>

#include <hdf5.h>
#include <mpi.h>

#include "types.hpp"

class BaseWriter {
public:
  BaseWriter(std::string filename, const MPI_Comm& comm = PETSC_COMM_WORLD);
  virtual ~BaseWriter() = default;

  virtual void write_image(const Image& image, const std::string& groupname) = 0;
  virtual void write_map(const Map& map, const std::string& groupname) = 0;

protected:
  MPI_Comm _comm;
  std::string filename;

  static const std::vector<std::string> _components;
};

#endif // BASEWRITER_HPP
