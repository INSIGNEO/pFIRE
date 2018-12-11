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

  intvector read_and_validate_image_header(const MPI_File &fh);
  intvector read_and_validate_mask_header(const MPI_File &fh);

  void copy_chunk_image(
      floating ***data, const std::vector<int> &subsize, const std::vector<int> &starts) const;
  void copy_chunk_mask(
      floating ***data, const std::vector<int> &subsize, const std::vector<int> &starts) const;
};

#endif // SHIRTLOADER_HPP
