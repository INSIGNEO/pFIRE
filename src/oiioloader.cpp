//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include "oiioloader.hpp"

#include "exceptions.hpp"
#include "file_utils.hpp"

//// ImageCache
// typedef and helpers for unique_ptr

struct ImageCacheDeleter {
  void operator()(OIIO::ImageCache *p) const
  {
    OIIO::ImageCache::destroy(p);
  }
};

ImageCache_unique create_unique_imagecache()
{
  return ImageCache_unique(OIIO::ImageCache::create());
}

const std::string OIIOLoader::loader_name = "OIIO";

ImageCache_unique OIIOLoader::cache(nullptr);

BaseLoader_unique OIIOLoader::Create_Loader(const std::string &path, MPI_Comm comm)
{
  return BaseLoader_unique(new OIIOLoader(path, comm));
}

OIIOLoader::OIIOLoader(const std::string &path, MPI_Comm comm) : BaseLoader(path, comm)
{
  if (!OIIOLoader::cache)
  {
    OIIOLoader::cache = create_unique_imagecache();
  }

  // open file and get metadata on image size
  OIIO::ImageSpec const *spec = OIIOLoader::cache->imagespec(OIIO::ustring(path));
  if (spec == nullptr)
  {
    throw_if_nonexistent(path);
    throw InvalidLoaderError(path);
  }
  this->_shape = {spec->width, spec->height, spec->depth};
}

void OIIOLoader::copy_scaled_chunk(
    floating ***data, const intvector &size, const intvector &corner_lo) const
{
  // petsc provides contiguous arrays so just get first element address
  floating *dataptr = &data[corner_lo[2]][corner_lo[1]][corner_lo[0]];
  intvector corner_hi(size.size(), 0);
  std::transform(size.cbegin(), size.cend(), corner_lo.cbegin(), corner_hi.begin(), std::plus<>());
  // OIIO scales floating point data between 0 and 1 internally on conversion
  bool res = cache->get_pixels(
      OIIO::ustring(_path), 0, 0, corner_lo[0], corner_hi[0], corner_lo[1], corner_hi[1],
      corner_lo[2], corner_hi[2], 0, 1, OIIO::TypeDesc::DOUBLE, dataptr, OIIO::AutoStride,
      OIIO::AutoStride, OIIO::AutoStride);
  if (!res)
  {
    throw_if_nonexistent(_path);
    throw InvalidLoaderError(_path);
  }
}
