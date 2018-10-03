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
  std::cout << "End of constructor" << std::endl;
}


void Elastic::autoregister()
{
  floatvector2d::const_reverse_iterator it = m_v_nodespacings.crbegin();
  while(it != m_v_nodespacings.rend())
  {
    std::cout << "Nodespacing " << (*it)[0] << std::endl;
    innerloop();
    std::advance(it, 1);
    m_v_nodespacings.erase(it.base());
    m_p_map = m_p_map->interpolate(*it);
  }
}

void Elastic::innerloop()
{
  // setup resolution specific solution storage (delta a, rvec)
  allocate_ephemeral_workspace();
  // calculate lambda for loop
  floating lambda = 1.0;
  for(integer inum=0; inum<m_max_iter; inum++)
  {
    std::cout << "Iteration " << inum << std::endl;
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
  perr = VecWAXPY(*m_vp_grads[0], -1.0, *m_fixed.local_vec(),
                  *m_p_registered->local_vec());CHKERRABORT(m_comm, perr);
  duplicate_single_grad_to_stacked(0);
  perr = MatMult(*m_p_tmat, *m_p_stacked_grads, *m_rhs_vec);

  // solve for delta a
  // update map
  // warp image
}

void Elastic::calculate_node_spacings(){
  const intvector& imshape = m_fixed.shape();
  floatvector currspc = m_v_final_nodespacing;
  m_v_nodespacings.push_back(currspc);
  while(all_true(currspc.begin(), currspc.end(), imshape.begin(), imshape.end(), std::less<>()))
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
    PetscErrorCode perr = VecDuplicate(*m_fixed.local_vec(), tmp_vec);CHKERRABORT(m_comm, perr);
    m_vp_grads.emplace_back(tmp_vec);
  }
  // create "global" vector for stack compatible with basis matrix
  // should be compatible with all map bases of this size
  PetscErrorCode perr = MatCreateVecs(*m_p_map->basis(), m_p_stacked_grads.get(),
                                      nullptr);CHKERRABORT(m_comm, perr);
}

void Elastic::allocate_ephemeral_workspace(){
  // allocate rhs vec and solution storage, use existing displacements in map
  PetscErrorCode perr;
  m_delta_vec = create_unique_vec();
  m_rhs_vec = create_unique_vec();
  perr = VecDuplicate(*m_p_map->m_displacements, m_delta_vec.get());CHKERRABORT(m_comm, perr);
  perr = VecDuplicate(*m_p_map->m_displacements, m_rhs_vec.get());CHKERRABORT(m_comm, perr);
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
  integer offset = 0;
  for(auto&& pp_grad : m_vp_grads)
  {
    //Create IS for target indices in stacked array, add to array for teardown
    IS* src_is = new IS;
    PetscErrorCode perr = ISCreateStride(m_comm, m_size, 0, 1, src_is);CHKERRABORT(m_comm, perr);
    m_vp_iss.emplace_back(src_is);

    IS* tgt_is = new IS;
    perr = ISCreateStride(m_comm, m_size, offset, 1, tgt_is);CHKERRABORT(m_comm, perr);
    m_vp_iss.emplace_back(tgt_is);
  
    std::cout << src_is << " : " << tgt_is << std::endl;
    // Create scatterer and add to array
    VecScatter* tmp_sctr = new VecScatter;
    perr = VecScatterCreate(*pp_grad, *src_is, *m_p_stacked_grads, *tgt_is, tmp_sctr);CHKERRABORT(m_comm, perr);
    m_vp_scatterers.emplace_back(tmp_sctr);

    offset += m_size;
  }
}

void Elastic::calculate_tmat()
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
    m_vp_grads[idim+1] = fd::gradient_to_global_unique(*(m_fixed.dmda()), *m_vp_grads[0], idim);
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
