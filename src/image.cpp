#include "image.hpp"

Image::Image(const intvector shape, MPI_Comm comm)
             : m_ndim(shape.size()), m_comm(comm), m_shape(shape) //const on shape causes copy assignment (c++11)
{
  if(m_shape.size() != 3)
  {
    if(m_shape.size() == 2)
    {
      m_shape.push_back(1);
    }
    else
    {
      throw std::runtime_error("image shape should be 2D or 3D");
    }
  }

  initialize_dmda();
  initialize_vectors();
}

Image::Image(const Image& image)
  : m_ndim(image.m_ndim), m_comm(image.m_comm), m_shape(image.m_shape),
    m_dmda(image.m_dmda)
{
  initialize_vectors();
}

void Image::initialize_dmda()
{
  integer dof_per_node = 1;
  integer stencil_width = 1;
  PetscErrorCode perr;

  // Make sure things get gracefully cleaned up
  m_dmda = std::shared_ptr<DM>(new DM, DMDestroy);
  perr = DMDACreate3d(m_comm,
                      DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, //BCs
                      DMDA_STENCIL_STAR, //stencil shape
                      m_shape[0], m_shape[1], m_shape[2], //global mesh shape
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, //ranks per dim
                      dof_per_node, stencil_width, //dof per node, stencil size
                      nullptr, nullptr, nullptr, //partition sizes nullptr -> petsc chooses
                      m_dmda.get());CHKERRABORT(m_comm, perr);

  perr = DMSetUp(*(m_dmda));CHKERRABORT(m_comm, perr);
}

void Image::initialize_vectors()
{
  PetscErrorCode perr;
  m_localvec = create_shared_vec();
  m_globalvec = create_shared_vec();
  perr = DMCreateGlobalVector(*m_dmda, m_globalvec.get());CHKERRABORT(m_comm, perr);
  perr = DMCreateLocalVector(*m_dmda, m_localvec.get());CHKERRABORT(m_comm, perr);
}

void Image::update_local_from_global()
{
  PetscErrorCode perr = DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);CHKERRABORT(PETSC_COMM_WORLD, perr);
}


Vec_unique Image::gradient(integer dim)
{
  // New global vec must be a duplicate of image global 
  Vec_shared grad = create_shared_vec();
  PetscErrorCode perr = VecDuplicate(*(m_globalvec), grad.get());CHKERRABORT(m_comm, perr);

  // Ensure we have up to date ghost cells
  DMGlobalToLocalBegin(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);
  // Nothing really useful to interleave
  DMGlobalToLocalEnd(*m_dmda, *m_globalvec, INSERT_VALUES, *m_localvec);

  return fd::gradient_to_global_unique(*m_dmda, *m_localvec, dim);
}



