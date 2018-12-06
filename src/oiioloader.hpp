#ifndef OIIOLOADER_HPP
#define OIIOLOADER_HPP

#include <OpenImageIO/imagecache.h>

#include "baseloader.hpp"
#include "types.hpp"

struct ImageCacheDeleter;
using ImageCache_unique = std::unique_ptr<OIIO::ImageCache, ImageCacheDeleter>;

class OIIOLoader: public BaseLoader {
public:
  static const std::string loader_name;

  OIIOLoader(const std::string &path, MPI_Comm comm = PETSC_COMM_WORLD);

  ~OIIOLoader() = default;

  void copy_scaled_chunk(floating ***data, const intvector &size, const intvector &offset) const;

  static BaseLoader_unique Create_Loader(const std::string &path, MPI_Comm comm);

private:
  static ImageCache_unique cache;
};

#endif // OIIOLOADER_HPP
