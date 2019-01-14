#ifndef HDFWRITER_HPP
#define HDFWRITER_HPP

#include <string>

#include <hdf5.h>
#include <mpi.h>

#include "types.hpp"
#include "basewriter.hpp"

class HDFWriter: public BaseWriter {
public:
  HDFWriter(std::string filename, const MPI_Comm& comm);
  ~HDFWriter();

  void write_image(const Image& image, const std::string& groupname);
  void write_map(const Map& map, const std::string& groupname);

  static const std::string writer_name;
  static const std::vector<std::string> extensions;

protected:
  std::string h5_filename;

private:
  hid_t _file_h;

  void create_or_truncate_h5();
};

#endif // HDFWRITER_HPP
