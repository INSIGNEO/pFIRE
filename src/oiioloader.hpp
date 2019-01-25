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
