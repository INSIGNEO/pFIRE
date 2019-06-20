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

#ifndef BASELOADER_HPP
#define BASELOADER_HPP

#include <functional>
#include <map>

#include "types.hpp"

class BaseLoader {
public:
  using loader_creator =
      std::function<BaseLoader_unique(const std::string& filename, MPI_Comm comm)>;

  using loader_map = std::map<std::string, loader_creator>;

  BaseLoader(const std::string& path, MPI_Comm comm = MPI_COMM_WORLD) : _comm(comm), _path(path) {}

  virtual ~BaseLoader() = default;

  const intcoord& shape() const { return _shape;}
  const std::string& path() const { return _path;}
  MPI_Comm comm() const { return _comm;}

  static const loader_map& loaders() { return *_loaders; }

  virtual void copy_scaled_chunk(floating*** data, const intvector& size,
                                 const intvector& offset) const = 0;

  static bool register_loader(const std::string& name, loader_creator loader);

  static BaseLoader_unique find_loader(const std::string& name, MPI_Comm comm = MPI_COMM_WORLD);

protected:
  void set_shape(const intcoord& shape){ _shape = shape;};

private:
  MPI_Comm _comm;
  std::string _path;
  intcoord _shape;

  static std::unique_ptr<loader_map> _loaders;
};

template <class T, class U>
void norm_convert(U* outptr, const T* inptr, uinteger size, double min, double max)
{
  U fscale = (U)(max - min);
  U fofs = (U)min;

  for (size_t idx = 0; idx < size; idx++)
  {
    outptr[idx] = (inptr[idx] - fofs) / fscale;
  }
}

#endif // BASELOADER_HPP
