#include "elastic.hpp"

Elastic::Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing)
  : m_comm(fixed.comm()), m_imgdims(fixed.ndim()), m_mapdims(m_imgdims+1), m_size(fixed.size()),
    m_iternum(0), m_fixed(fixed), m_moved(moved), m_v_final_nodespacing(nodespacing),
    m_p_stacked_grads(create_unique_vec())
{
  // TODO: image compatibility checks (maybe write Image.iscompat(Image foo)
  
  // make sure nodespacing is compatible with image
  if(m_fixed.ndim() != m_v_final_nodespacing.size())
  {
    throw std::runtime_error("number of nodespacings must match number of image dimensions");
  }

  // work out intermediate node spacings
  calculate_node_spacings();

  // make map, need to ensure basis is always the same layout 
  m_p_map = std::make_unique<Map>(fixed, m_v_nodespacings.back());

  // not in initializer to avoid copy until we know images are compatible and we can proceed
  m_p_registered = moved.copy();

  // set scratchpad storage, scatterers:
  allocate_persistent_workspace();
  create_scatterers();
}


void Elastic::autoregister()
{
  integer loop_cnt = 1;
  PetscPrintf(m_comm, "Beginning elastic registration\n");
  if(m_imgdims == 2)
  {
    PetscPrintf(m_comm, "Target nodespacing: %.1f, %.1f\n",
                m_v_final_nodespacing[0], m_v_final_nodespacing[1]);
  }
  else
  {
    PetscPrintf(m_comm, "Target node spacing: %.1f, %.1f, %.1f\n",
                m_v_final_nodespacing[0], m_v_final_nodespacing[1], m_v_final_nodespacing[2]);
  }
  PetscPrintf(m_comm, "Using %i generations\n\n", m_v_nodespacings.size());

  auto it = m_v_nodespacings.crbegin();
  while(it != m_v_nodespacings.rend())
  {
    PetscPrintf(m_comm, "Generation %i, ", loop_cnt++);
    if(m_imgdims == 2)
    {
      PetscPrintf(m_comm, "Nodespacing %.1f, %.1f\n", (*it)[0], (*it)[1]);
    } 
    else
    {
      PetscPrintf(m_comm, "Nodespacing %.1f, %.1f, %1.f\n", (*it)[0], (*it)[1], (*it)[2]);
    }

    innerloop();
    std::advance(it, 1);
    if(it == m_v_nodespacings.rend())
    {
      break;
    }
    m_v_nodespacings.erase(it.base());
    m_p_map = m_p_map->interpolate(*it);
  }
}

void Elastic::innerloop()
{
  // setup map resolution specific solution storage (tmat, delta a, rvec)
  allocate_ephemeral_workspace();
  // calculate lambda for loop
  floating lambda = 1.0;
  for(integer inum=0; inum<m_max_iter; inum++)
  {
    PetscPrintf(m_comm, "Iteration %i:\n", inum);
    innerstep(lambda);
    //check convergence and break
  }
}

void Elastic::innerstep(floating lambda)
{
  // calculate up to date tmat
  calculate_tmat();

  // calculate tmat2 and precondition
  Mat_unique normmat = create_unique_mat();
  // TODO: can we reuse here?
  PetscErrorCode perr = MatTransposeMatMult(*m_p_tmat, *m_p_tmat, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
                                            normmat.get());CHKERRABORT(m_comm, perr);
  // precondition tmat2
  // TODO

  // calculate tmat2 + lambda*lapl2
  perr = MatAXPY(*normmat, lambda, *m_p_map->laplacian(),
                 DIFFERENT_NONZERO_PATTERN);CHKERRABORT(m_comm, perr);

  // calculate rvec, to do this need to reuse stacked vector for [f-m f-m f-m f-m]
  perr = VecWAXPY(*m_vp_grads[0], -1.0, *m_fixed.global_vec(),
                  *m_p_registered->global_vec());CHKERRABORT(m_comm, perr);
  duplicate_single_grad_to_stacked(0);
  perr = MatMultTranspose(*m_p_tmat, *m_p_stacked_grads, *m_rhs_vec);
  
  // solve for delta a
  KSP_unique m_ksp = create_unique_ksp();
  KSPCreate(m_comm, m_ksp.get());
  KSPSetOperators(*m_ksp, *normmat, *normmat);
  KSPSolve(*m_ksp, *m_rhs_vec, *m_delta_vec); 
  // update map
  m_p_map->update(*m_delta_vec);
  // warp image
//  m_p_registered = m_p_map->warp(m_moved);
}

void Elastic::calculate_node_spacings()
{
  const intvector& imshape = m_fixed.shape();
  floatvector currspc = m_v_final_nodespacing;
  m_v_nodespacings.push_back(currspc);
  while(all_true_varlen(currspc.begin(), currspc.end(), imshape.begin(), imshape.end(), std::less<>()))
  {
    std::transform(currspc.begin(), currspc.end(), currspc.begin(),
                   [](floating a) -> floating{return a*2;});
    m_v_nodespacings.push_back(currspc);
  }
}

void Elastic::allocate_persistent_workspace()
{
  //create "local" vectors for gradient storage, one per map dim
  for(integer idim=0; idim < m_mapdims; idim++)
  {
    Vec* tmp_vec = new Vec;
    PetscErrorCode perr = VecDuplicate(*m_fixed.global_vec(), tmp_vec);CHKERRABORT(m_comm, perr);
    m_vp_grads.emplace_back(tmp_vec);
    integer gsize, lsize, ostart, oend;
    VecGetSize(*tmp_vec, &gsize);
    VecGetLocalSize(*tmp_vec, &lsize);
    VecGetOwnershipRange(*tmp_vec, &ostart, &oend);
  }

  // create local vector to allow gradient calculation
  m_gradlocal = create_unique_vec();
  PetscErrorCode perr = VecDuplicate(*m_fixed.local_vec(), m_gradlocal.get());CHKERRABORT(m_comm, perr);

  // create "global" vector for stack compatible with basis matrix
  // should be compatible with all map bases of this size
  perr = MatCreateVecs(*m_p_map->basis(), nullptr,
                                      m_p_stacked_grads.get());CHKERRABORT(m_comm, perr);
  integer gsize, lsize, ostart, oend;
  VecGetSize(*m_p_stacked_grads, &gsize);
  VecGetLocalSize(*m_p_stacked_grads, &lsize);
  VecGetOwnershipRange(*m_p_stacked_grads, &ostart, &oend);
}

void Elastic::allocate_ephemeral_workspace(){
  // allocate rhs vec and solution storage, use existing displacements in map
  PetscErrorCode perr;
  m_delta_vec = create_unique_vec();
  m_rhs_vec = create_unique_vec();
  perr = VecDuplicate(*m_p_map->m_displacements, m_delta_vec.get());CHKERRABORT(m_comm, perr);
  perr = VecDuplicate(*m_p_map->m_displacements, m_rhs_vec.get());CHKERRABORT(m_comm, perr);

  // setup map resolution specific tmat storage, can use basis layout
  m_p_tmat = create_unique_mat();
  perr = MatDuplicate(*m_p_map->m_basis, MAT_DO_NOT_COPY_VALUES, m_p_tmat.get());CHKERRABORT(m_comm, perr);

}


void Elastic::scatter_grads_to_stacked()
{
  for(integer idim=0; idim < m_mapdims; idim++)
  {
    PetscErrorCode perr = VecScatterBegin(*m_vp_scatterers[idim], *m_vp_grads[idim], *m_p_stacked_grads,
                                          INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(*m_vp_scatterers[idim], *m_vp_grads[idim], *m_p_stacked_grads,
                                        INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
  }
}

void Elastic::duplicate_single_grad_to_stacked(size_t idx)
{
  for(integer idim=0; idim < m_mapdims; idim++)
  {
    PetscErrorCode perr = VecScatterBegin(*m_vp_scatterers[idim], *m_vp_grads[idx], *m_p_stacked_grads,
                                          INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
    perr = VecScatterEnd(*m_vp_scatterers[idim], *m_vp_grads[idx], *m_p_stacked_grads,
                                        INSERT_VALUES, SCATTER_FORWARD);CHKERRABORT(m_comm, perr);
  }
}


void Elastic::create_scatterers()
{
  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed obj.
  PetscErrorCode perr = DMDAGetAO(*m_fixed.dmda(), &ao_petsctonat); // Destroying this would break the dmda

  integer offset = 0; // Track offset
  for(auto&& pp_grad : m_vp_grads)
  {
    // Get extents of local data in grad array
    integer startelem, datasize;
    perr = VecGetOwnershipRange(*pp_grad, &startelem, &datasize);CHKERRABORT(m_comm, perr);
    datasize -= startelem;

    // Create IS for source indices in stacked array, keep in vector for memory management
    m_vp_iss.push_back(create_unique_is());
    IS* p_src_is = m_vp_iss.back().get(); //grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem, 1, p_src_is);CHKERRABORT(m_comm, perr);
    // Convert with AO to map from petsc to natural ordering (do here because tgt_is is offset)
    perr = AOApplicationToPetscIS(ao_petsctonat, *m_vp_iss.back());CHKERRABORT(m_comm, perr);

    m_vp_iss.push_back(create_unique_is());
    IS* p_tgt_is = m_vp_iss.back().get(); //grab a temp copy to avoid vector accesses below
    perr = ISCreateStride(m_comm, datasize, startelem+offset, 1, p_tgt_is);CHKERRABORT(m_comm, perr);

    // Create scatterer and add to array
    m_vp_scatterers.push_back(create_unique_vecscatter());
    perr = VecScatterCreate(*pp_grad, *p_src_is, *m_p_stacked_grads, *p_tgt_is, m_vp_scatterers.back().get());CHKERRABORT(m_comm, perr);

    offset += m_size;
  }
}

void Elastic::calculate_tmat()
{
  // Calculate average intensity 0.5(f+m)
  // Constant offset needed later, does not affect gradients
  PetscErrorCode perr = VecSet(*m_vp_grads[0], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);
  //NB Z = aX + bY + cZ has call signature VecAXPBYPCZ(Z, a, b, c, X, Y) because reasons....
  perr = VecAXPBYPCZ(*m_vp_grads[0], 0.5, 0.5, 1,  *m_fixed.global_vec(),
                     *m_moved.global_vec());CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter this to local for later
  perr = DMGlobalToLocalBegin(*m_fixed.dmda(), *m_vp_grads[0], INSERT_VALUES, *m_gradlocal);
  perr = DMGlobalToLocalEnd(*m_fixed.dmda(), *m_vp_grads[0], INSERT_VALUES, *m_gradlocal);

  // find average gradients
  for(integer idim=0; idim < m_fixed.ndim(); idim++)
  {
    //likely change this to avoid mallocs/frees
    fd::gradient_existing(*(m_fixed.dmda()), *m_gradlocal, *m_vp_grads[idim+1], idim);
  }

  // Negate average intensity to get 1 - 0.5(f+m) as needed by algorithm
  perr = VecScale(*m_vp_grads[0], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter grads into stacked vector
  scatter_grads_to_stacked();

  // 3. copy basis into p_tmat, can do trivially because we already duplicated in setup
  perr = MatCopy(*m_p_map->basis(), *m_p_tmat, SAME_NONZERO_PATTERN);CHKERRABORT(m_comm, perr);

  // 4. left diagonal multiply p_tmat with stacked vector
  perr = MatDiagonalScale(*m_p_tmat, *m_p_stacked_grads, nullptr);CHKERRABORT(m_comm, perr);
}
