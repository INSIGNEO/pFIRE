#include "elastic.hpp"

#include "fd_routines.hpp"
#include "iterator_routines.hpp"

Elastic::Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing)
  : m_comm(fixed.comm()), m_imgdims(fixed.ndim()), m_mapdims(m_imgdims+1), m_size(fixed.size()),
    m_iternum(0), m_fixed(fixed), m_moved(moved), m_v_final_nodespacing(nodespacing)
{
  // TODO: image compatibility checks (maybe write Image.iscompat(Image foo)
  // TODO: enforce normalization
  
  // make sure nodespacing is compatible with image
  if(m_fixed.ndim() != m_v_final_nodespacing.size())
  {
    throw std::runtime_error("number of nodespacings must match number of image dimensions");
  }

  #ifdef VERBOSEDEBUG
  PetscPrintf(m_comm, "Fixed image initial contents:\n");
  VecView(*m_fixed.global_vec(), PETSC_VIEWER_STDOUT_(m_comm));
  PetscPrintf(m_comm, "Moved image initial contents:\n");
  VecView(*m_fixed.global_vec(), PETSC_VIEWER_STDOUT_(m_comm));
  #endif

  // not in initializer to avoid copy until we know images are compatible and we can proceed
  m_p_registered = moved.copy();

  // work out intermediate node spacings
  calculate_node_spacings();

  // make map, need to ensure basis is always the same layout 
  m_p_map = std::make_unique<Map>(fixed, m_v_nodespacings.back());

  // set scratchpad storage, scatterers:
  m_workspace = std::make_shared<WorkSpace>(fixed, *m_p_map);


}


void Elastic::autoregister()
{
  integer loop_count = 1;
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
    PetscPrintf(m_comm, "Generation %i, ", loop_count);
    if(m_imgdims == 2)
    {
      PetscPrintf(m_comm, "Nodespacing %.1f, %.1f\n", (*it)[0], (*it)[1]);
    } 
    else
    {
      PetscPrintf(m_comm, "Nodespacing %.1f, %.1f, %1.f\n", (*it)[0], (*it)[1], (*it)[2]);
    }
    innerloop(loop_count);
    std::advance(it, 1);
    if(it == m_v_nodespacings.rend())
    {
      break;
    }
    m_v_nodespacings.erase(it.base());
    m_p_map = m_p_map->interpolate(*it);
    m_workspace->reallocate_ephemeral_workspace(*m_p_map);
    m_p_registered = m_p_map->warp(m_moved, *m_workspace);
    loop_count++;
  }
}

void Elastic::innerloop(integer outer_count)
{
  // setup map resolution specific solution storage (tmat, delta a, rvec)
  // calculate lambda for loop
  #ifdef DEBUG_DUMP_INTERMEDIATES
  save_debug_frame(outer_count, 0);
  #endif //DEBUG_DUMP_INTERMEDIATES

  floating lambda = 20.0;
  for(integer inum=1; inum<=m_max_iter; inum++)
  {
    PetscPrintf(m_comm, "Iteration %i:\n", inum);
    innerstep(lambda);
   
    #ifdef DEBUG_DUMP_INTERMEDIATES
    save_debug_frame(outer_count, inum);
    #endif //DEBUG_DUMP_INTERMEDIATES

    //check convergence and break if below threshold
    floating posmax, negmax;
    PetscErrorCode perr = VecMax(*m_workspace->m_delta, nullptr, &posmax);CHKERRABORT(m_comm, perr);
    perr = VecMin(*m_workspace->m_delta, nullptr, &negmax);CHKERRABORT(m_comm, perr);
    floating amax = std::max(std::fabs(posmax), std::fabs(negmax));
    PetscPrintf(m_comm, "Maximum displacement: %.2f\n", amax);
    if(amax < m_convergence_thres)
    {
      break;
    }
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
  perr = VecWAXPY(*m_workspace->m_globaltmps[0], -1.0, *m_p_registered->global_vec(),
                 *m_fixed.global_vec());CHKERRABORT(m_comm, perr);
  m_workspace->duplicate_single_grad_to_stacked(0);
  perr = MatMultTranspose(*m_workspace->m_tmat, *m_workspace->m_stacktmp, *m_workspace->m_rhs);
  
  block_precondition();
  // solve for delta a
  KSP_unique m_ksp = create_unique_ksp();
  KSPCreate(m_comm, m_ksp.get());
  KSPSetOperators(*m_ksp, *normmat, *normmat);
  KSPSetFromOptions(*m_ksp);
  KSPSetUp(*m_ksp);
  KSPSolve(*m_ksp, *m_workspace->m_rhs, *m_workspace->m_delta); 
  // update map
  m_p_map->update(*m_workspace->m_delta);
  // warp image
  m_p_registered = m_p_map->warp(m_moved, *m_workspace);
  m_p_registered->normalize();
}

void Elastic::calculate_node_spacings()
{
  const intvector& imshape = m_fixed.shape();
  floatvector currspc = m_v_final_nodespacing;
  m_v_nodespacings.push_back(currspc);
  while(all_true_varlen(currspc.begin(), currspc.end(), imshape.begin(), imshape.end(),
        [](floating x, integer y) -> bool {return (y/x) > 2.0;}))
  {
    std::transform(currspc.begin(), currspc.end(), currspc.begin(),
                   [](floating a) -> floating{return a*2;});
    m_v_nodespacings.push_back(currspc);
  }

  PetscPrintf(m_comm, "Calculated spacings:\n");
  for(auto spcs = m_v_nodespacings.crbegin(); spcs != m_v_nodespacings.crend(); spcs++)
  {
    PetscPrintf(m_comm, "Spacing: ");
    for(auto spc = spcs->cbegin(); spc != spcs->cend(); spc++)
    {
      PetscPrintf(m_comm, "%g ", *spc);
    }
    PetscPrintf(m_comm, "\n");
  }
}

void Elastic::calculate_tmat()
{
  // Calculate average intensity 0.5(f+m)
  // Constant offset needed later, does not affect gradients
  PetscErrorCode perr = VecSet(*m_workspace->m_globaltmps[m_fixed.ndim()], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);
  //NB Z = aX + bY + cZ has call signature VecAXPBYPCZ(Z, a, b, c, X, Y) because reasons....
  perr = VecAXPBYPCZ(*m_workspace->m_globaltmps[m_fixed.ndim()], 0.5, 0.5, 1,  *m_fixed.global_vec(),
                     *m_p_registered->global_vec());CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter this to local for later
  perr = DMGlobalToLocalBegin(*m_fixed.dmda(), *m_workspace->m_globaltmps[m_fixed.ndim()], INSERT_VALUES,
                              *m_workspace->m_localtmp);CHKERRABORT(m_comm, perr);
  perr = DMGlobalToLocalEnd(*m_fixed.dmda(), *m_workspace->m_globaltmps[m_fixed.ndim()], INSERT_VALUES,
                            *m_workspace->m_localtmp);CHKERRABORT(m_comm, perr);

  // find average gradients
  for(integer idim=0; idim < m_fixed.ndim(); idim++)
  {
    //likely change this to avoid mallocs/frees
    fd::gradient_existing(*(m_fixed.dmda()), *m_workspace->m_localtmp,
                          *m_workspace->m_globaltmps[idim], idim);
  }

  // Negate average intensity to get 1 - 0.5(f+m) as needed by algorithm
  perr = VecScale(*m_workspace->m_globaltmps[m_fixed.ndim()], -1.0);CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter grads into stacked vector
  m_workspace->scatter_grads_to_stacked();

  // 3. copy basis into p_tmat, can do trivially because we already duplicated in setup
  perr = MatCopy(*m_p_map->basis(), *m_workspace->m_tmat, SAME_NONZERO_PATTERN);CHKERRABORT(m_comm, perr);

  // 4. left diagonal multiply p_tmat with stacked vector
  perr = MatDiagonalScale(*m_workspace->m_tmat, *m_workspace->m_stacktmp, nullptr);CHKERRABORT(m_comm, perr);
}

void Elastic::block_precondition()
{
  // Normalize luminance block of matrix to spatial blocks using diagonal norm
  Vec_unique diag = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(*normmat, diag.get(), nullptr);CHKERRABORT(m_comm, perr);
  perr = MatGetDiagonal(*normmat, *diag);CHKERRABORT(m_comm, perr);
  
  // index where spatial dims stop and intensity dim starts
  integer crit_idx = m_p_map->size() * m_p_map->m_ndim;

  // Find rank-local sums for spatial and luminance blocks
  integer rowstart, rowend;
  perr = VecGetOwnershipRange(*diag, &rowstart, &rowend);CHKERRABORT(m_comm, perr);
  integer localsize = rowend - rowstart;
  floating norm[2]; // norm[0] is spatial, norm[1] is luminance
  integer spt_end = crit_idx - rowstart;
  spt_end = (spt_end < 0) ? 0 : (spt_end > localsize) ? localsize : spt_end;

  floating *ptr;
  perr = VecGetArray(*diag, &ptr);CHKERRABORT(m_comm, perr);
  for(integer idx=0; idx < spt_end; idx++)
  {
    norm[0] += ptr[idx];
  }
  for(integer idx=spt_end; idx < localsize; idx++)
  {
    norm[1] += ptr[idx];
  }
  perr = VecRestoreArray(*diag, &ptr);CHKERRABORT(m_comm, perr);
  
  // MPI_AllReduce to sum over all processes
  MPI_Allreduce(MPI_IN_PLACE, norm, 2, MPI_DOUBLE, MPI_SUM, m_comm); 

  // calculate average of norms and scaling factor
  norm[0] /= crit_idx;
  norm[1] /= m_p_map->size();
  floating lum_scale = norm[0]/norm[1];

  // reuse diag vector to hold scaling values
  perr = VecGetArray(*diag, &ptr);CHKERRABORT(m_comm, perr);
  for(integer idx=0; idx < spt_end; idx++)
  {
    ptr[idx] = 1;
  }
  for(integer idx=spt_end; idx < localsize; idx++)
  {
    ptr[idx] = lum_scale;
  }
  perr = VecRestoreArray(*diag, &ptr);CHKERRABORT(m_comm, perr);

  perr = MatDiagonalScale(*normmat, *diag, nullptr);CHKERRABORT(m_comm, perr);
}

void Elastic::save_debug_frame(integer ocount, integer icount)
{
  std::ostringstream fname;
  fname << "dbg"
        << std::setfill('0') << std::setw(2) << ocount << "_" 
        << std::setw(3) << icount << ".png";
  m_p_registered->save_OIIO(fname.str());
}
