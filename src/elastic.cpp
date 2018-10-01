#include "elastic.hpp"

Elastic::Elastic(const Image& fixed, const Image& moved, const floatvector node_spacing)
  : m_comm(fixed.comm()), m_imgdims(fixed.ndim()), m_mapdims(m_imgdims+1),
    m_size(fixed.size()), m_fixed(fixed), m_moved(moved)
{
  // TODO: image compatibility checks (maybe write Image.iscompat(Image foo)
  
  // make map, need to ensure basis is always the same layout 

  // required scratchpad storage:
  //  tmat -> duplicate of map.basis()
  //  stacked vector compatible with tmat
  //  local work vectors compatible with image
  //  scatters from locals to stacked above
  //  map must have laplacian, basis and a vectors compatible
}

void Elastic::create_scratch()
{
  //create "local" vectors for gradient storage, one per map dim
  for(integer idim=0; idim < m_mapdims; idim++)
  {
    Vec* tmp_vec = new Vec;
    PetscErrorCode perr = VecDuplicate(*m_fixed.local_vec(), tmp_vec);CHKERRABORT(m_comm, perr);
    m_vp_iss.emplace_back(tmp_vec);
  }
  // create "global" vector for stack compatible with basis matrix
  // should be compatible with all map bases of this size
  PetscErrorCode perr = MatCreateVecs(*m_p_map->basis(), m_p_stacked_grads.get(),
                                      nullptr);CHKERRABORT(m_comm, perr);
}

void Elastic::do_scatter_to_stacked()
{
  for(integer idim=0; idim < m_mapdims; idim++)
  {
    PetscErrorCode perr = VecScatterBegin(*m_vp_scatterers[idim], *m_vp_grads[idim], *m_p_stacked_grads,
                                          INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
  }
  for(integer idim=0; idim < m_mapdims; idim++)
  {
    PetscErrorCode perr = VecScatterEnd(*m_vp_scatterers[idim], *m_vp_grads[idim], *m_p_stacked_grads,
                                        INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
  }
}

void Elastic::create_scatterers()
{
  integer offset = 0;
  for(auto&& pp_grad : m_vp_grads)
  {
    //Create IS for target indices in stacked array, add to array for teardown
    IS* tmp_is = new IS;
    PetscErrorCode perr = ISCreateStride(m_comm, m_size, offset, 1, tmp_is);CHKERRABORT(m_comm, perr);
    m_vp_iss.emplace_back(tmp_is);

    // Create scatterer and add to array
    VecScatter* tmp_sctr = new VecScatter;
    perr = VecScatterCreate(*pp_grad, nullptr, *m_p_stacked_grads, *tmp_is, tmp_sctr);CHKERRABORT(m_comm, perr);
    m_vp_scatterers.emplace_back(tmp_sctr);

    offset += m_size;
  }
}

void Elastic::create_t_matrix()
{
  // Calculate average intensity 0.5(f+m)
  // Constant offset needed later, does not affect gradients
  PetscErrorCode perr = VecSet(*(m_vp_grads[0]), -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);
  //NB Z = aX + bY + cZ has call signature VecAXPBYPCZ(Z, a, b, c, X, Y) because reasons....
  perr = VecAXPBYPCZ(*(m_vp_grads[0]), 0.5, 0.5, 1,  *(m_fixed.local_vec()),
                     *(m_moved.local_vec()));CHKERRABORT(PETSC_COMM_WORLD, perr);

  // find average gradients
  for(integer idim=0; idim < m_fixed.ndim(); idim++)
  {
    //likely change this to avoid mallocs/frees
    m_vp_grads[idim+1] = fd::gradient_to_global_unique(*(fixed.dmda()), *vp_grads[0], idim);
  }

  // Negate average intensity to get 1 - 0.5(f+m) as needed by algorithm
  perr = VecScale(*vp_grads[0], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);

  // 2. scatter grads into stacked vector
  // 3. copy basis into p_tmat 
  // 4. diagonal multiply p_tmat with stacked vector

  return p_tmat;
}
