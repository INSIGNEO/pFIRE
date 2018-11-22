#include "oiioloader.hpp"

#include<OpenImageIO/imageio.h>

//// ImageCache
// typedef and helpers for unique_ptr
struct ImageCacheDeleter{void operator()(OIIO::ImageCache* p) const{OIIO::ImageCache::destroy(p);}};
using ImageCache_unique = std::unique_ptr<OIIO::ImageCache, ImageCacheDeleter>;
inline ImageCache_unique create_unique_imagecache()
{
  return ImageCache_unique(OIIO::ImageCache::create());
}


ImageCache_unique OIIOLoader::cache(nullptr);

const std::string OIIOLoader::loader_name = "OIIO";

BaseLoader_unique OIIOLoader::Create_Loader(const std::string &path, MPI_Comm comm)
{
  return BaseLoader_unique(new OIIOLoader(path, comm));
}


OIIOLoader::OIIOLoader(const std::string& path, MPI_Comm comm)
  : BaseLoader(path, comm)
{
  if(!OIIOLoader::cache)
  {
    OIIOLoader::cache = create_unique_imagecache();
  }

  // open file and get metadata on image size
  OIIO::ImageSpec const* spec = OIIOLoader::cache->imagespec(OIIO::ustring(path));
  if(spec == nullptr){
    throw std::runtime_error("Failed to open image file");
  }
  this->_shape = {spec->width, spec->height, spec->depth};
}

void OIIOLoader::copy_scaled_chunk(floating ***data, const intvector &size, 
                                  const intvector &corner_lo) const
{
  // petsc provides contiguous arrays so just get first element address
  floating *dataptr = &data[corner_lo[2]][corner_lo[1]][corner_lo[0]];
  intvector corner_hi(size.size(), 0);
  std::transform(size.cbegin(), size.cend(), corner_lo.cbegin(), corner_hi.begin(), std::plus<>());
  // OIIO scales floating point data between 0 and 1 internally on conversion
  bool res = cache->get_pixels(OIIO::ustring(_path), 0, 0, 
                               corner_lo[0], corner_hi[0],
                               corner_lo[1], corner_hi[1],
                               corner_lo[2], corner_hi[2],
                               0, 1,
                               OIIO::TypeDesc::DOUBLE, dataptr,
                               OIIO::AutoStride, OIIO::AutoStride, OIIO::AutoStride);
  if(!res)
  {
    throw std::runtime_error("Failed to load image");
  }
}
