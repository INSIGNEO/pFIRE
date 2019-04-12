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

#include "workspace.hpp"

WorkSpace::WorkSpace(const Image& image, const Map& map)
    : m_comm(image.comm()), m_dmda(image.dmda()), m_size(image.size()),
      m_globaltmps(std::vector<Vec_unique>()), m_iss(std::vector<IS_unique>()),
      m_r_scatterers(std::vector<VecScatter_unique>()),
      m_nr_scatterers(std::vector<VecScatter_unique>()), m_stacktmp(create_unique_vec()),
      m_localtmp(create_unique_vec()), m_delta(create_unique_vec()), m_rhs(create_unique_vec()),
      m_tmat(create_unique_mat()), ephemeral_count(0)
{
  // create "local" vectors for gradient storage, one per map dim
  for (uinteger idim = 0; idim < image.ndim() + 1; idim++)
  {
    Vec_unique tmp_vec = create_unique_vec();
    PetscErrorCode perr = VecDuplicate(*image.global_vec(), tmp_vec.get());
    debug_creation(*tmp_vec, std::string("grads_") + std::to_string(static_cast<int>(idim)));
    CHKERRABORT(m_comm, perr);
    m_globaltmps.push_back(std::move(tmp_vec));
  }

  // create local vector to allow gradient calculation
  m_localtmp = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*image.local_vec(), m_localtmp.get());
  CHKERRABORT(m_comm, perr);

  // create "global" vector for stack compatible with basis matrix
  // should be compatible with all map bases of this size
  m_stacktmp = create_unique_vec();
  perr = MatCreateVecs(*map.basis(), nullptr, m_stacktmp.get());
  debug_creation(*m_stacktmp, "workspace vector");
  CHKERRABORT(m_comm, perr);

  create_reordering_scatterers();
  create_nonreordering_scatterers();
  reallocate_ephemeral_workspace(map);
}

void WorkSpace::reallocate_ephemeral_workspace(const Map& map)
{
  ephemeral_count++;
  // allocate rhs vec and solution storage, use existing displacements in map
  PetscErrorCode perr;
  m_delta = create_unique_vec();
  m_rhs = create_unique_vec();
  perr = VecDuplicate(*map.m_displacements, m_delta.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*m_delta, std::string("solution_storage_") + std::to_string(ephemeral_count));
  perr = VecSet(*m_delta, 0.);
  CHKERRABORT(m_comm, perr);
  perr = VecDuplicate(*map.m_displacements, m_rhs.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*m_delta, std::string("rhs_storage_") + std::to_string(ephemeral_count));
  perr = VecSet(*m_rhs, 0.);
  CHKERRABORT(m_comm, perr);
}

void WorkSpace::scatter_stacked_to_grads()
{
  for (size_t idim = 0; idim < m_r_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(
        *m_r_scatterers[idim], *m_stacktmp, *m_globaltmps[idim], INSERT_VALUES, SCATTER_REVERSE);
    CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(
        *m_r_scatterers[idim], *m_stacktmp, *m_globaltmps[idim], INSERT_VALUES, SCATTER_REVERSE);
    CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::scatter_grads_to_stacked()
{
  for (size_t idim = 0; idim < m_r_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(
        *m_r_scatterers[idim], *m_globaltmps[idim], *m_stacktmp, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(
        *m_r_scatterers[idim], *m_globaltmps[idim], *m_stacktmp, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::duplicate_single_grad_to_stacked(size_t idx)
{
  for (size_t idim = 0; idim < m_r_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(
        *m_r_scatterers[idim], *m_globaltmps[idx], *m_stacktmp, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(
        *m_r_scatterers[idim], *m_globaltmps[idx], *m_stacktmp, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::scatter_stacked_to_grads_noreorder()
{
  for (size_t idim = 0; idim < m_nr_scatterers.size(); idim++)
  {
    PetscErrorCode perr = VecScatterBegin(
        *m_nr_scatterers[idim], *m_stacktmp, *m_globaltmps[idim], INSERT_VALUES, SCATTER_REVERSE);
    CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(
        *m_nr_scatterers[idim], *m_stacktmp, *m_globaltmps[idim], INSERT_VALUES, SCATTER_REVERSE);
    CHKERRABORT(m_comm, perr);
  }
}

void WorkSpace::create_nonreordering_scatterers()
{
  integer offset = 0; // Track offset
  for (auto&& pp_grad : m_globaltmps)
  {
    // Get extents of local data in grad array
    integer startelem, datasize;
    PetscErrorCode perr = VecGetOwnershipRange(*pp_grad, &startelem, &datasize);
    CHKERRABORT(m_comm, perr);
    datasize -= startelem;

    // Create IS for source indices in stacked array, keep in vector for memory management
    m_iss.push_back(create_unique_is());
    IS* p_src_is = m_iss.back().get(); // grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem, 1, p_src_is);
    CHKERRABORT(m_comm, perr);

    m_iss.push_back(create_unique_is());
    IS* p_tgt_is = m_iss.back().get(); // grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem + offset, 1, p_tgt_is);
    CHKERRABORT(m_comm, perr);

    // Create scatterer and add to array
    m_nr_scatterers.push_back(create_unique_vecscatter());
    perr =
        VecScatterCreate(*pp_grad, *p_src_is, *m_stacktmp, *p_tgt_is, m_nr_scatterers.back().get());
    CHKERRABORT(m_comm, perr);

    offset += m_size;
  }
}

void WorkSpace::create_reordering_scatterers()
{
  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed
                    // obj.
  PetscErrorCode perr = DMDAGetAO(*m_dmda, &ao_petsctonat); // Destroying this would break the dmda
  CHKERRABORT(m_comm, perr);

  integer offset = 0; // Track offset
  for (auto&& pp_grad : m_globaltmps)
  {
    // Get extents of local data in grad array
    integer startelem, datasize;
    perr = VecGetOwnershipRange(*pp_grad, &startelem, &datasize);
    CHKERRABORT(m_comm, perr);
    datasize -= startelem;

    // Create IS for source indices in stacked array, keep in vector for memory management
    m_iss.push_back(create_unique_is());
    IS* p_src_is = m_iss.back().get(); // grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem, 1, p_src_is);
    CHKERRABORT(m_comm, perr);
    // Convert with AO to map from petsc to natural ordering (do here because tgt_is is offset)
    perr = AOApplicationToPetscIS(ao_petsctonat, *m_iss.back());
    CHKERRABORT(m_comm, perr);

    m_iss.push_back(create_unique_is());
    IS* p_tgt_is = m_iss.back().get(); // grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem + offset, 1, p_tgt_is);
    CHKERRABORT(m_comm, perr);

    // Create scatterer and add to array
    m_r_scatterers.push_back(create_unique_vecscatter());
    perr =
        VecScatterCreate(*pp_grad, *p_src_is, *m_stacktmp, *p_tgt_is, m_r_scatterers.back().get());
    CHKERRABORT(m_comm, perr);

    offset += m_size;
  }
}
