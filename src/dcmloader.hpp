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
