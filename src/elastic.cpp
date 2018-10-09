#include "elastic.hpp"

Elastic::Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing)
  : m_comm(fixed.comm()), m_imgdims(fixed.ndim()), m_mapdims(m_imgdims+1), m_size(fixed.size()),
    m_iternum(0), m_fixed(fixed), m_moved(moved), m_v_final_nodespacing(nodespacing)
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
  m_workspace = std::make_shared<WorkSpace>(fixed, *m_p_map);
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
  m_workspace->reallocate_ephemeral_workspace(*m_p_map);
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
  PetscErrorCode perr = MatTransposeMatMult(*m_workspace->m_tmat, *m_workspace->m_tmat,
                                            MAT_INITIAL_MATRIX, PETSC_DEFAULT,
                                            normmat.get());CHKERRABORT(m_comm, perr);
  // precondition tmat2
  // TODO

  // calculate tmat2 + lambda*lapl2
  perr = MatAXPY(*normmat, lambda, *m_p_map->laplacian(),
                 DIFFERENT_NONZERO_PATTERN);CHKERRABORT(m_comm, perr);

  // calculate rvec, to do this need to reuse stacked vector for [f-m f-m f-m f-m]
  perr = VecWAXPY(*m_workspace->m_globaltmps[0], -1.0, *m_fixed.global_vec(),
                  *m_p_registered->global_vec());CHKERRABORT(m_comm, perr);
  m_workspace->duplicate_single_grad_to_stacked(0);
  perr = MatMultTranspose(*m_workspace->m_tmat, *m_workspace->m_stacktmp, *m_workspace->m_rhs);
  
  // solve for delta a
  KSP_unique m_ksp = create_unique_ksp();
  KSPCreate(m_comm, m_ksp.get());
  KSPSetOperators(*m_ksp, *normmat, *normmat);
  KSPSolve(*m_ksp, *m_workspace->m_rhs, *m_workspace->m_delta); 
  // update map
  m_p_map->update(*m_workspace->m_delta);
  // warp image
  //m_p_registered = m_p_map->warp(m_moved, *m_workspace->m_stacktmp);
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

void Elastic::calculate_tmat()
{
  // Calculate average intensity 0.5(f+m)
  // Constant offset needed later, does not affect gradients
  PetscErrorCode perr = VecSet(*m_workspace->m_globaltmps[0], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);
  //NB Z = aX + bY + cZ has call signature VecAXPBYPCZ(Z, a, b, c, X, Y) because reasons....
  perr = VecAXPBYPCZ(*m_workspace->m_globaltmps[0], 0.5, 0.5, 1,  *m_fixed.global_vec(),
                     *m_moved.global_vec());CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter this to local for later
  perr = DMGlobalToLocalBegin(*m_fixed.dmda(), *m_workspace->m_globaltmps[0], INSERT_VALUES,
                              *m_workspace->m_localtmp);CHKERRABORT(m_comm, perr);
  perr = DMGlobalToLocalEnd(*m_fixed.dmda(), *m_workspace->m_globaltmps[0], INSERT_VALUES,
                            *m_workspace->m_localtmp);CHKERRABORT(m_comm, perr);

  // find average gradients
  for(integer idim=0; idim < m_fixed.ndim(); idim++)
  {
    //likely change this to avoid mallocs/frees
    fd::gradient_existing(*(m_fixed.dmda()), *m_workspace->m_localtmp,
                          *m_workspace->m_globaltmps[idim+1], idim);
  }

  // Negate average intensity to get 1 - 0.5(f+m) as needed by algorithm
  perr = VecScale(*m_workspace->m_globaltmps[0], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter grads into stacked vector
  m_workspace->scatter_grads_to_stacked();

  // 3. copy basis into p_tmat, can do trivially because we already duplicated in setup
  perr = MatCopy(*m_p_map->basis(), *m_workspace->m_tmat, SAME_NONZERO_PATTERN);CHKERRABORT(m_comm, perr);

  // 4. left diagonal multiply p_tmat with stacked vector
  perr = MatDiagonalScale(*m_workspace->m_tmat, *m_workspace->m_stacktmp, nullptr);CHKERRABORT(m_comm, perr);
}
