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

#include "image.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <petscdmda.h>
#include <petscvec.h>

#include "fd_routines.hpp"
#include "iterator_routines.hpp"
#include "map.hpp"
#include "indexing.hpp"
#include "math_utils.hpp"

#include "baseloader.hpp"

// Public Methods

Image::Image(const intvector& shape, MPI_Comm comm)
  : ImageBase(shape, comm)
{
}

std::unique_ptr<Image> Image::duplicate(const ImageBase& img)
{
  return std::unique_ptr<Image>(new Image(img));
}

std::unique_ptr<Image> Image::copy(const ImageBase& img)
{
  // private copy c'tor prohibits use of std::make_unique, would otherwise do:
  // std::unique_ptr<Image> new_img = std::make_unique<Image>(*this);
  std::unique_ptr<Image> new_img(new Image(img));

  new_img->copy_data(img);
  return new_img;
}

std::unique_ptr<Image>
Image::load_file(const std::string& path, const ImageBase* existing, MPI_Comm comm)
{
  BaseLoader_unique loader = BaseLoader::find_loader(path, comm);

  // if image passed assert sizes match and duplicate, otherwise create new image given size
  std::unique_ptr<Image> new_image;
  if (existing != nullptr)
  {
    comm = existing->comm();
    if (!all_true(
            loader->shape().begin(), loader->shape().end(), existing->shape().begin(),
            existing->shape().end(), std::equal_to<>()))
    {
      throw std::runtime_error("New image must have same shape as existing");
    }
    new_image = Image::duplicate(*existing);
  }
  else
  {
    new_image = std::make_unique<Image>(loader->shape(), comm);
  }

  intvector shape(3, 0), offset(3, 0);
  PetscErrorCode perr = DMDAGetCorners(
      *new_image->dmda(), &offset[0], &offset[1], &offset[2], &shape[0], &shape[1], &shape[2]);
  CHKERRABORT(comm, perr);
  // std::transform(shape.cbegin(), shape.cend(), offset.cbegin(), shape.begin(), std::minus<>());

  floating*** vecptr(nullptr);
  perr = DMDAVecGetArray(*new_image->dmda(), *new_image->global_vec(), &vecptr);
  CHKERRABORT(comm, perr);
  loader->copy_scaled_chunk(vecptr, shape, offset);
  perr = DMDAVecRestoreArray(*new_image->dmda(), *new_image->global_vec(), &vecptr);
  CHKERRABORT(comm, perr);

  new_image->normalize();

  return new_image;
}

floating Image::masked_normalize(const Mask& mask)
{
  Vec_unique tmp = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*m_globalvec, tmp.get());CHKERRABORT(m_comm, perr);
  perr = VecPointwiseMult(*tmp, *mask.global_vec(), *m_globalvec);CHKERRABORT(m_comm, perr);

  floating norm;
  perr = VecSum(*tmp, &norm);CHKERRABORT(m_comm, perr);
  norm = mask.npoints() / norm;
  perr = VecScale(*m_globalvec, norm);CHKERRABORT(m_comm, perr);

  return norm;
}

// Protected Methods

Image::Image(const ImageBase& image)
  : ImageBase(image)
{
}
