#ifndef HDFWRITER_HPP
#define HDFWRITER_HPP

#include <string>

#include <hdf5.h>
#include <mpi.h>

#include "types.hpp"

class HDFWriter {
public:
  HDFWriter(const std::string& filename, const MPI_Comm& comm, bool truncate_existing = true);
  ~HDFWriter();

  void write_image(const Image& image, const std::string& groupname);

  void write_map(const Map& map, const std::string& groupname);

private:
  hid_t create_or_open(std::string filename, hid_t file_props);
  hid_t create_or_truncate(std::string filename, hid_t file_props);

  MPI_Comm _comm;
  std::string _filename;
  hid_t _file_h;

  static const std::vector<std::string> _components;
};

#endif
