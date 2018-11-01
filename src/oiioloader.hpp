#ifndef OIIOLOADER_HPP
#define OIIOLOADER_HPP

#include "types.hpp"
#include "baseloader.hpp"

class OIIOLoader: public BaseLoader{
public:

  static const std::string loader_name;

  OIIOLoader(const std::string& path, MPI_Comm comm=PETSC_COMM_WORLD);

  ~OIIOLoader() = default;

  void copy_scaled_chunk(floating ***data, const intvector &size, const intvector &offset) const;

  static BaseLoader_unique Create_Loader(const std::string& path, MPI_Comm comm);

private:

  static ImageCache_unique cache;
};

#endif //OIIOLOADER_HPP
