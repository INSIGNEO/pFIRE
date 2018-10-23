#include "workspace.hpp"

WorkSpace::WorkSpace(const Image& image, const Map& map)
  : m_comm(image.comm()), m_dmda(image.dmda()), m_size(image.size()),
    m_globaltmps(std::vector<Vec_unique>()), m_iss(std::vector<IS_unique>()),
    m_scatterers(std::vector<VecScatter_unique>()), m_stacktmp(create_unique_vec()),
    m_localtmp(create_unique_vec()), m_delta(create_unique_vec()), m_rhs(create_unique_vec()),
    m_tmat(create_unique_mat())
{
  //create "local" vectors for gradient storage, one per map dim
  for(uinteger idim=0; idim < image.ndim()+1; idim++)
  {
    Vec_unique tmp_vec = create_unique_vec();
    PetscErrorCode perr = VecDuplicate(*image.global_vec(), tmp_vec.get());CHKERRABORT(m_comm, perr);
    m_globaltmps.push_back(std::move(tmp_vec));
  }

  // create local vector to allow gradient calculation
  m_localtmp = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*image.local_vec(), m_localtmp.get());CHKERRABORT(m_comm, perr);

  // create "global" vector for stack compatible with basis matrix
  // should be compatible with all map bases of this size
  m_stacktmp = create_unique_vec();
  perr = MatCreateVecs(*map.basis(), nullptr, m_stacktmp.get());CHKERRABORT(m_comm, perr);

  create_scatterers();
  reallocate_ephemeral_workspace(map);
}

void WorkSpace::reallocate_ephemeral_workspace(const Map& map){
  // allocate rhs vec and solution storage, use existing displacements in map
  PetscErrorCode perr;
  m_delta = create_unique_vec();
  m_rhs = create_unique_vec();
  perr = VecDuplicate(*map.m_displacements, m_delta.get());CHKERRABORT(m_comm, perr);
  perr = VecDuplicate(*map.m_displacements, m_rhs.get());CHKERRABORT(m_comm, perr);

  // setup map resolution specific tmat storage, can use basis layout
  m_tmat = create_unique_mat();
  perr = MatDuplicate(*map.m_basis, MAT_DO_NOT_COPY_VALUES, m_tmat.get());CHKERRABORT(m_comm, perr);

}

void WorkSpace::scatter_stacked_to_grads()
{
  for(size_t idim=0; idim < m_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(*m_scatterers[idim], *m_stacktmp, *m_globaltmps[idim],
                                          INSERT_VALUES, SCATTER_REVERSE);CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(*m_scatterers[idim], *m_stacktmp, *m_globaltmps[idim],
                                        INSERT_VALUES, SCATTER_REVERSE);CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::scatter_grads_to_stacked()
{
  for(size_t idim=0; idim < m_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(*m_scatterers[idim], *m_globaltmps[idim], *m_stacktmp,
                                          INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(*m_scatterers[idim], *m_globaltmps[idim], *m_stacktmp,
                                        INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::duplicate_single_grad_to_stacked(size_t idx)
{
  for(size_t idim=0; idim < m_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(*m_scatterers[idim], *m_globaltmps[idx], *m_stacktmp,
                                          INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(*m_scatterers[idim], *m_globaltmps[idx], *m_stacktmp,
                                        INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::create_scatterers()
{
  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed obj.
  PetscErrorCode perr = DMDAGetAO(*m_dmda, &ao_petsctonat); // Destroying this would break the dmda

  integer offset = 0; // Track offset
  for(auto&& pp_grad : m_globaltmps)
  {
    // Get extents of local data in grad array
    integer startelem, datasize;
    perr = VecGetOwnershipRange(*pp_grad, &startelem, &datasize);CHKERRABORT(m_comm, perr);
    datasize -= startelem;

    // Create IS for source indices in stacked array, keep in vector for memory management
    m_iss.push_back(create_unique_is());
    IS* p_src_is = m_iss.back().get(); //grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem, 1, p_src_is);CHKERRABORT(m_comm, perr);
    // Convert with AO to map from petsc to natural ordering (do here because tgt_is is offset)
    perr = AOApplicationToPetscIS(ao_petsctonat, *m_iss.back());CHKERRABORT(m_comm, perr);

    m_iss.push_back(create_unique_is());
    IS* p_tgt_is = m_iss.back().get(); //grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem+offset, 1, p_tgt_is);CHKERRABORT(m_comm, perr);

    // Create scatterer and add to array
    m_scatterers.push_back(create_unique_vecscatter());
    perr = VecScatterCreate(*pp_grad, *p_src_is, *m_stacktmp, *p_tgt_is,
                            m_scatterers.back().get());CHKERRABORT(m_comm, perr);

    offset += m_size;
  }
}
