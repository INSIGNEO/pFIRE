#ifndef DCMLOADER_HPP
#define DCMLOADER_HPP

#include <dcmtk/config/osconfig.h>
#include <dcmtk/dcmdata/dctk.h>

#include "baseloader.hpp"
#include "types.hpp"

class DCMLoader: public BaseLoader {
public:
  static const std::string loader_name;

  DCMLoader(const std::string &path, MPI_Comm comm = PETSC_COMM_WORLD);

  ~DCMLoader() = default;

  void copy_scaled_chunk(floating ***data, const intvector &size, const intvector &offset) const;

  static BaseLoader_unique Create_Loader(const std::string &path, MPI_Comm comm);

private:
  mutable DcmFileFormat _datafile;
};

#endif // DCMLOADER_HPP
