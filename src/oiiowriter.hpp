#ifndef OIIOWRITER_HPP
#define OIIOWRITER_HPP

#include <string>

#include "basewriter.hpp"

class OIIOWriter : public BaseWriter {
  public:
  OIIOWriter(std::string filename, const MPI_Comm& comm = PETSC_COMM_WORLD);
  ~OIIOWriter() = default;

  void write_image(const Image& image);
  
  void write_map(const Map& map);
};

#endif // OIIOWRITER_HPP
