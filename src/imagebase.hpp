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

#ifndef IMAGEBASE_HPP
#define IMAGEBASE_HPP

#include <petscdmda.h>

#include "types.hpp"

class ImageBase {
public:
  explicit ImageBase(const intvector& shape, MPI_Comm comm = PETSC_COMM_WORLD);

  static std::unique_ptr<Image> duplicate(const ImageBase& img);
  static std::unique_ptr<Image> copy(const ImageBase& img);

  // Allow RO access to member variables
  // Note that datavec and dmda remain mutable in this way
  MPI_Comm comm() const { return m_comm; }
  uinteger ndim() const { return m_ndim; }
  const intvector& shape() const { return m_shape; }
  integer size() const { return m_shape[0] * m_shape[1] * m_shape[2]; }
  std::shared_ptr<const Vec> global_vec() const { return m_globalvec; }
  std::shared_ptr<const Vec> local_vec() const { return m_localvec; }
  std::shared_ptr<DM> dmda() const { return m_dmda; }

  const floating* get_raw_data_ro() const;
  void release_raw_data_ro(const floating*& ptr) const;
  void copy_data(const ImageBase& img);
  Vec_unique get_raw_data_row_major() const;
  Vec_unique get_raw_data_natural() const;

  Vec_unique gradient(integer dim);
  floating normalize();
  floating masked_normalize(const Mask& mask);
  void binarize();

  template <typename inttype>
  std::vector<inttype> mpi_get_offset() const;
  template <typename inttype>
  std::vector<inttype> mpi_get_chunksize() const;

  Vec_unique scatter_to_zero() const;

/*!
   * Calculate mutual information between this and a second image.
   *
   * The two images must share a DMDA.
   *
   * @param other the second image for mutual information calculation
   * @return The mutual information between images
   */
  floating mutual_information(const ImageBase &other);

protected:
  void update_local_from_global();

  explicit ImageBase(const ImageBase& image);
  ImageBase& operator=(const ImageBase& image);

  MPI_Comm m_comm;
  uinteger m_ndim;
  intvector m_shape;
  Vec_shared m_localvec, m_globalvec;
  DM_shared m_dmda;
  //  std::shared_ptr<Mask> mask;

  void initialize_dmda();
  void initialize_vectors();

  integer instance_id;
  static integer instance_id_counter;

private:
  const size_t mi_resolution = 100;
};

template <typename inttype>
std::vector<inttype> ImageBase::mpi_get_chunksize() const
{
  // This routine only makes sense to use for integer types
  static_assert(std::is_integral<inttype>::value, "Integral element type required");

  intvector sizes(3, 0);
  PetscErrorCode perr =
      DMDAGetCorners(*m_dmda, nullptr, nullptr, nullptr, &sizes[0], &sizes[1], &sizes[2]);
  CHKERRXX(perr);

  std::vector<inttype> out(sizes.begin(), sizes.end());

  return out;
}

template <typename inttype>
std::vector<inttype> ImageBase::mpi_get_offset() const
{
  // This routine only makes sense to use for integer types
  static_assert(std::is_integral<inttype>::value, "Integral element type required");

  intvector offsets(3, 0);
  PetscErrorCode perr =
      DMDAGetCorners(*m_dmda, &offsets[0], &offsets[1], &offsets[2], nullptr, nullptr, nullptr);
  CHKERRXX(perr);

  std::vector<inttype> out(offsets.begin(), offsets.end());

  return out;
}

#endif // IMAGEBASE_HPP
