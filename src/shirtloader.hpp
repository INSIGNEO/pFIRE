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

#ifndef SHIRTLOADER_HPP
#define SHIRTLOADER_HPP

#include "baseloader.hpp"
#include "types.hpp"

class ShIRTLoader: public BaseLoader {
public:
  static const std::string loader_name;

  ShIRTLoader(const std::string &path, MPI_Comm comm = PETSC_COMM_WORLD);

  ~ShIRTLoader() = default;

  void copy_scaled_chunk(floating ***data, const intvector &size, const intvector &offset) const;

  static BaseLoader_unique Create_Loader(const std::string &path, MPI_Comm comm);

private:
  enum shirt_file_type : integer { undef = 0, image = 1, mask = 2 };

  shirt_file_type _file_type;

  using shirt_header_dtype = int32_t;

  using image_data_dtype = float;
  static constexpr MPI_Datatype image_data_mpi_type = MPI_FLOAT;
  static constexpr integer image_header_length = 4;
  static constexpr integer image_header_bytes = image_header_length * sizeof(shirt_header_dtype);

  using mask_data_dtype = int16_t;
  static constexpr MPI_Datatype mask_data_mpi_type = MPI_SHORT;
  static constexpr integer mask_header_length = 3;
  static constexpr integer mask_header_bytes = mask_header_length * sizeof(shirt_header_dtype);

  intcoord read_and_validate_image_header(const MPI_File &fh);
  intcoord read_and_validate_mask_header(const MPI_File &fh);

  void copy_chunk_image(
      floating ***data, const std::vector<int> &subsize, const std::vector<int> &starts) const;
  void copy_chunk_mask(
      floating ***data, const std::vector<int> &subsize, const std::vector<int> &starts) const;
};

#endif // SHIRTLOADER_HPP
