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

#include "mask.hpp"

#include "baseloader.hpp"
#include "iterator_routines.hpp"
#include "exceptions.hpp"

Mask::Mask(const intvector& shape, MPI_Comm comm)
  : ImageBase(shape, comm)
{
}

std::unique_ptr<Mask> Mask::duplicate(const ImageBase& img)
{
  return std::unique_ptr<Mask>(new Mask(img));
}

std::unique_ptr<Mask> Mask::copy(const ImageBase& img)
{
  // private copy c'tor prohibits use of std::make_unique, would otherwise do:
  // std::unique_ptr<Mask> new_img = std::make_unique<Mask>(*this);
  std::unique_ptr<Mask> new_img(new Mask(img));

  new_img->copy_data(img);
  return new_img;
}

std::unique_ptr<Mask> Mask::full_image(const ImageBase& img)
{
  std::unique_ptr<Mask> new_mask = Mask::duplicate(img);
  PetscErrorCode perr = VecSet(*new_mask->m_globalvec, 1.0);
  CHKERRXX(perr);
  
  return new_mask;
}

std::unique_ptr<Mask>
Mask::load_file(const std::string& path, const ImageBase* existing, MPI_Comm comm)
{
  BaseLoader_unique loader = BaseLoader::find_loader(path, comm);

  // if mask passed assert sizes match and duplicate, otherwise create new mask given size
  std::unique_ptr<Mask> new_mask;
  if (existing != nullptr)
  {
    comm = existing->comm();
    if (!all_true(
            loader->shape().begin(), loader->shape().end(), existing->shape().begin(),
            existing->shape().end(), std::equal_to<>()))
    {
      throw InternalError("New mask must have same shape as existing", __FILE__, __LINE__);
    }
    new_mask = Mask::duplicate(*existing);
  }
  else
  {
    new_mask = std::make_unique<Mask>(loader->shape(), comm);
  }

  intvector shape(3, 0), offset(3, 0);
  PetscErrorCode perr = DMDAGetCorners(
      *new_mask->dmda(), &offset[0], &offset[1], &offset[2], &shape[0], &shape[1], &shape[2]);
  CHKERRXX(perr);
  // std::transform(shape.cbegin(), shape.cend(), offset.cbegin(), shape.begin(), std::minus<>());

  floating*** vecptr(nullptr);
  perr = DMDAVecGetArray(*new_mask->dmda(), *new_mask->global_vec(), &vecptr);
  CHKERRXX(perr);
  loader->copy_scaled_chunk(vecptr, shape, offset);
  perr = DMDAVecRestoreArray(*new_mask->dmda(), *new_mask->global_vec(), &vecptr);
  CHKERRXX(perr);

  new_mask->binarize();

  return new_mask;
}

integer Mask::npoints() const
{
  floating total_points;
  PetscErrorCode perr = VecSum(*m_globalvec, &total_points);
  CHKERRXX(perr);

  return static_cast<integer>(total_points);
}

// Protected Methods

Mask::Mask(const ImageBase& mask)
  : ImageBase(mask)
{
}
