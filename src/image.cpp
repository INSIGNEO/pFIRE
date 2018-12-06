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

#ifdef USE_OIIO
#include <OpenImageIO/imageio.h>
#endif // USE_OIIO

#include "baseloader.hpp"

// Public Methods

Image::Image(const intvector& shape, MPI_Comm comm)
    : m_comm(comm), m_ndim(shape.size()),
      m_shape(shape), // const on shape causes copy assignment (c++11)
      m_localvec(create_unique_vec()), m_globalvec(create_unique_vec()), m_dmda(create_shared_dm())
{
  if (m_shape.size() != 3)
  {
    if (m_shape.size() == 2)
    {
      m_shape.push_back(1);
    }
    else
    {
      throw std::runtime_error("image shape should be 2D or 3D");
    }
  }
  if (m_shape[2] == 1)
  {
    m_ndim = 2;
  }

  initialize_dmda();
  initialize_vectors();
}

std::unique_ptr<Image> Image::duplicate() const
{
  return std::unique_ptr<Image>(new Image(*this));
}

std::unique_ptr<Image> Image::copy() const
{
  // private copy c'tor prohibits use of std::make_unique, would otherwise do:
  // std::unique_ptr<Image> new_img = std::make_unique<Image>(*this);
  std::unique_ptr<Image> new_img(new Image(*this));

  PetscErrorCode perr = VecCopy(*m_localvec, *new_img->m_localvec);
  CHKERRABORT(m_comm, perr);
  perr = VecCopy(*m_globalvec, *new_img->m_globalvec);
  CHKERRABORT(m_comm, perr);

  return new_img;
}

std::unique_ptr<Image>
Image::load_file(const std::string& path, const Image* existing, MPI_Comm comm)
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
    new_image = existing->duplicate();
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

  return new_image;
}
// Return scale factor
floating Image::normalize()
{
  floating norm;
  PetscErrorCode perr = VecSum(*m_globalvec, &norm);
  CHKERRABORT(m_comm, perr);
  norm = this->size() / norm;
  perr = VecScale(*m_globalvec, norm);
  CHKERRABORT(m_comm, perr);
  return norm;
}
/*
void Image::set_mask(std::shared_ptr<Mask> mask)
{
  //TODO iscompat()
  this->mask = mask;
}

void Image::masked_normalize()
{
  if(m_mask.get() == nullptr)
  {
    throw std::runtime_error("mask must be set first")
  }

  Vec_unique tmp = create_unique_vec();
  PetscErrorcode perr = VecDuplicate(*m_globalvec, tmp.get());CHKERRABORT(m_comm, perr);
  perr = VecPointwiseMult(*tmp, *this->mask->global_vec(), *m_globalvec);CHKERRABORT(m_comm, perr);

  floating norm;
  perr = VecSum(*tmp, &norm);CHKERRABORT(m_comm, perr);
  norm = this->mask->npoints() / norm;
  perr = VecScale(*m_globalvec, norm);CHKERRABORT(m_comm, perr);

  return norm;
}*/

#ifdef USE_OIIO
void Image::save_OIIO(std::string filename)
{
  Vec_unique imgvec = create_unique_vec();
  PetscErrorCode perr = DMDACreateNaturalVector(*m_dmda, imgvec.get());
  CHKERRABORT(m_comm, perr);
  perr = DMDAGlobalToNaturalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *imgvec);
  CHKERRABORT(m_comm, perr);
  perr = DMDAGlobalToNaturalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *imgvec);
  CHKERRABORT(m_comm, perr);

  if (m_comm != MPI_COMM_SELF)
  {
    imgvec = scatter_to_zero(*imgvec);
  }

  integer rank;
  MPI_Comm_rank(m_comm, &rank);
  if (rank == 0)
  {
    OIIO::ImageOutput* img = OIIO::ImageOutput::create(filename);
    if (img == nullptr)
    {
      throw std::runtime_error("Failed to open image output file");
    }
    OIIO::ImageSpec spec(m_shape[0], m_shape[1], 1, OIIO::TypeDesc::UINT16);
    img->open(filename, spec);

    floating* pixdata;
    perr = VecGetArray(*imgvec, &pixdata);
    CHKERRABORT(m_comm, perr);
    img->write_image(OIIO::TypeDesc::DOUBLE, pixdata);
    img->close();
    perr = VecRestoreArray(*imgvec, &pixdata);
    CHKERRABORT(m_comm, perr);

    OIIO::ImageOutput::destroy(img);
  }
  MPI_Barrier(m_comm);
}
#endif // USE_OIIO

// Protected Methods

Image::Image(const Image& image)
    : m_comm(image.m_comm), m_ndim(image.m_ndim), m_shape(image.m_shape),
      m_localvec(create_shared_vec()), m_globalvec(create_shared_vec()), m_dmda(image.m_dmda)
{
  initialize_vectors();
}

void Image::initialize_dmda()
{
  integer dof_per_node = 1;
  integer stencil_width = 1;
  PetscErrorCode perr;

  // Make sure things get gracefully cleaned up
  m_dmda = create_shared_dm();
  perr = DMDACreate3d(
      m_comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, // BCs
      DMDA_STENCIL_STAR,                                                     // stencil shape
      m_shape[0], m_shape[1], m_shape[2],                                    // global mesh shape
      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,                              // ranks per dim
      dof_per_node, stencil_width, // dof per node, stencil size
      nullptr, nullptr, nullptr,   // partition sizes nullptr -> petsc chooses
      m_dmda.get());
  CHKERRABORT(m_comm, perr);

  perr = DMSetUp(*(m_dmda));
  CHKERRABORT(m_comm, perr);
}

void Image::initialize_vectors()
{
  PetscErrorCode perr;
  m_localvec = create_shared_vec();
  m_globalvec = create_shared_vec();
  perr = DMCreateGlobalVector(*m_dmda, m_globalvec.get());
  CHKERRABORT(m_comm, perr);
  perr = DMCreateLocalVector(*m_dmda, m_localvec.get());
  CHKERRABORT(m_comm, perr);
}

void Image::update_local_from_global()
{
  PetscErrorCode perr = DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
}

Vec_unique Image::gradient(integer dim)
{
  // New global vec must be a duplicate of image global
  Vec_shared grad = create_shared_vec();
  PetscErrorCode perr = VecDuplicate(*(m_globalvec), grad.get());
  CHKERRABORT(m_comm, perr);

  // Ensure we have up to date ghost cells
  DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  // Nothing really useful to interleave
  DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);

  return fd::gradient_to_global_unique(*m_dmda, *m_localvec, dim);
}

Vec_unique Image::scatter_to_zero(Vec& vec)
{
  Vec_unique new_vec = create_unique_vec();
  VecScatter_unique sct = create_unique_vecscatter();
  VecScatterCreateToZero(vec, sct.get(), new_vec.get());
  VecScatterBegin(*sct, vec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(*sct, vec, *new_vec, INSERT_VALUES, SCATTER_FORWARD);

  return new_vec;
}

const floating* Image::get_raw_data_ro() const
{
  const floating* ptr;
  PetscErrorCode perr = VecGetArrayRead(*m_globalvec, &ptr);
  CHKERRABORT(m_comm, perr);
  return ptr;
}

void Image::release_raw_data_ro(const floating*& ptr) const
{
  PetscErrorCode perr = VecRestoreArrayRead(*m_globalvec, &ptr);
  CHKERRABORT(m_comm, perr);
}
