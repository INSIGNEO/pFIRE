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

#include "types.hpp"

#include "exceptions.hpp"
#include "gridvariable.hpp"
#include "baseloader.hpp"
#include "iterator_routines.hpp"

class ImageBase: public GridVariable {
public:
  ImageBase(const intcoord& shape, const MPI_Comm& comm);

  template <typename ImageClass>
  static std::shared_ptr<
      typename std::enable_if<std::is_base_of<ImageBase, ImageClass>::value, ImageClass>::type>
  load_image(std::string image_path, const ImageClass *existing = nullptr,
             MPI_Comm comm = MPI_COMM_WORLD);

  floating normalize();

  Vec_unique difference_global(const ImageBase& other) const;

  const Vec& data_vector() const { return this->global_vector();}

  floating mutual_information(const ImageBase &other);

protected:
};

template <typename ImageClass>
std::shared_ptr<typename std::enable_if<std::is_base_of<ImageBase, ImageClass>::value, ImageClass>::type>
ImageBase::load_image(std::string image_path, const ImageClass *existing, MPI_Comm comm)
{
  BaseLoader_unique loader = BaseLoader::find_loader(image_path, comm);

  // if image passed assert sizes match and duplicate, otherwise create new image given size
  std::shared_ptr<ImageClass> new_image;
  if (existing != nullptr)
  {
    comm = existing->comm();
    if (!all_true(loader->shape().begin(), loader->shape().end(), existing->shape().begin(),
                  existing->shape().end(), std::equal_to<>()))
    {
      throw InternalError("New image must have same shape as existing", __FILE__, __LINE__);
    }
    new_image = GridVariable::duplicate(*existing);
  }
  else
  {
    new_image = std::make_shared<ImageClass>(loader->shape(), comm);
  }

  intvector shape(3, 0), offset(3, 0);
  PetscErrorCode perr = DMDAGetCorners(new_image->dmda(), &offset[0], &offset[1], &offset[2],
                                       &shape[0], &shape[1], &shape[2]);
  CHKERRXX(perr);
  // std::transform(shape.cbegin(), shape.cend(), offset.cbegin(), shape.begin(), std::minus<>());

  floating ***vecptr(nullptr);
  perr = DMDAVecGetArray(new_image->dmda(), new_image->global_vector(), &vecptr);
  CHKERRXX(perr);
  loader->copy_scaled_chunk(vecptr, shape, offset);
  perr = DMDAVecRestoreArray(new_image->dmda(), new_image->global_vector(), &vecptr);
  CHKERRXX(perr);

  auto extents = new_image->minmax();
  PetscPrintf(comm, "Loader minmax: %f, %f\n", extents.first, extents.second);

  new_image->normalize();
  new_image->update_local_vector();

  return new_image;
}

#endif // IMAGEBASE_HPP
