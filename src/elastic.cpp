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

#include "elastic.hpp"

#include <iomanip>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "basewriter.hpp"
#include "fd_routines.hpp"
#include "file_utils.hpp"
#include "infix_iterator.hpp"
#include "iterator_routines.hpp"
#include "math_utils.hpp"
#include "petsc_debug.hpp"

namespace bf = boost::filesystem;

Elastic::Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing,
    const ConfigurationBase& configuration)
  : m_comm(fixed.comm()), configuration(configuration), m_imgdims(fixed.ndim()),
    m_mapdims(m_imgdims + 1), m_size(fixed.size()), m_iternum(0), m_fixed(fixed), m_moved(moved),
    m_v_nodespacings(floatvector2d()), m_v_final_nodespacing(nodespacing),
    m_p_registered(std::shared_ptr<Image>(nullptr)), m_p_map(std::unique_ptr<Map>(nullptr)),
    m_workspace(std::shared_ptr<WorkSpace>(nullptr)), normmat(create_unique_mat())
{
  // TODO: image compatibility checks (maybe write Image.iscompat(Image foo)
  // TODO: enforce normalization

  // make sure nodespacing is compatible with image
  if (m_fixed.ndim() != m_v_final_nodespacing.size())
  {
    throw std::runtime_error("number of nodespacings must match number of image dimensions");
  }
  if (m_v_final_nodespacing.size() == 2)
  {
    m_v_final_nodespacing.push_back(1);
  }

  m_max_iter = configuration.grab<integer>("max_iterations"); 

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
  std::ostringstream nsmsg;
  nsmsg << "Target nodespacing: ";
  std::copy_n(
      m_v_final_nodespacing.cbegin(), m_imgdims, infix_ostream_iterator<integer>(nsmsg, " "));
  nsmsg << std::endl;
  PetscPrintf(m_comm, nsmsg.str().c_str());
  PetscPrintf(m_comm, "Using %i generations\n\n", m_v_nodespacings.size());

  auto it = m_v_nodespacings.crbegin();
  while (it != m_v_nodespacings.rend())
  {
    nsmsg << "Nodespacing: ";
    std::copy_n(it->cbegin(), m_imgdims, infix_ostream_iterator<integer>(nsmsg, " "));
    nsmsg << std::endl;
    PetscPrintf(m_comm, nsmsg.str().c_str());

    innerloop(loop_count);
    std::advance(it, 1);
    if (it == m_v_nodespacings.rend())
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
  if (configuration.grab<bool>("save_intermediate_frames"))
  {
    save_debug_frame(outer_count, 0);
  }

  bool recalculate_lambda = false;
  try
  {
    m_lambda = configuration.grab<floating>("lambda");
  }
  catch (std::invalid_argument& err)
  {
    m_lambda = Elastic::k_lambda_default;
    recalculate_lambda = true;
  }
  for (integer inum = 1; inum <= m_max_iter; inum++)
  {
    PetscPrintf(m_comm, "Iteration %i:\n", inum);
    innerstep(inum, recalculate_lambda);
    recalculate_lambda = false;

    if (configuration.grab<bool>("save_intermediate_frames"))
    {
      save_debug_frame(outer_count, inum);
    }

    // check convergence and break if below threshold
    floating posmax, negmax;
    PetscErrorCode perr = VecMax(*m_workspace->m_delta, nullptr, &posmax);
    CHKERRABORT(m_comm, perr);
    perr = VecMin(*m_workspace->m_delta, nullptr, &negmax);
    CHKERRABORT(m_comm, perr);
    floating amax = std::max(std::fabs(posmax), std::fabs(negmax));
    PetscPrintf(m_comm, "Maximum displacement: %.2f\n", amax);
    floating aavg;
    perr = VecNorm(*m_workspace->m_delta, NORM_2, &aavg);
    aavg /= m_p_map->size();
    PetscPrintf(m_comm, "Average displacement: %.2f\n", aavg);
    if (aavg < m_convergence_thres)
    {
      PetscPrintf(m_comm, "Generation %i converged after %i iterations.\n\n", outer_count, inum);
      break;
    }
  }
}

void Elastic::innerstep(integer inum, bool recalculate_lambda)
{
  floating lambda_mult = configuration.grab<floating>("lambda_mult");
  // calculate up to date tmat
  calculate_tmat(inum);

  // calculate tmat2 and precondition
  normmat = create_unique_mat();
  // TODO: can we reuse here?
  PetscErrorCode perr = MatTransposeMatMult(*m_workspace->m_tmat, *m_workspace->m_tmat,
      MAT_INITIAL_MATRIX, PETSC_DEFAULT, normmat.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*normmat, std::string("Mat_normal") + std::to_string(inum));
  // precondition tmat2
  block_precondition(*normmat, m_p_map->size(), m_p_map->m_ndim);

  floating lapldiagavg = diagonal_sum(*m_p_map->laplacian()) / (m_p_map->size() * m_p_map->m_ndim);
  floating tdiagavg = diagonal_sum(*normmat) / m_p_map->size();
  floating lapl_mult = tdiagavg / lapldiagavg;
  PetscPrintf(m_comm, "laplavg = %g, tavg = %g, mult=%g\n", lapldiagavg, tdiagavg, lapl_mult);

  if (recalculate_lambda)
  {
    m_lambda = approximate_optimum_lambda(
        *normmat, *m_p_map->laplacian(), lapl_mult, m_lambda, 10.0, 30, k_lambda_min);
    PetscPrintf(m_comm, "Calculated smoothing factor: %.2f\n", m_lambda);
  }

  floating total_mult = lapl_mult * lambda_mult * m_lambda;
  // calculate tmat2 + lambda*lapl2
  perr = MatAXPY(*normmat, total_mult, *m_p_map->laplacian(), DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(m_comm, perr);

  // calculate rvec, to do this need to reuse stacked vector for [f-m f-m f-m f-m]
  perr = VecWAXPY(
      *m_workspace->m_globaltmps[0], -1.0, *m_p_registered->global_vec(), *m_fixed.global_vec());
  CHKERRABORT(m_comm, perr);
  m_workspace->duplicate_single_grad_to_stacked(0);
  perr = MatMultTranspose(*m_workspace->m_tmat, *m_workspace->m_stacktmp, *m_workspace->m_rhs);
  CHKERRABORT(m_comm, perr);

  // apply memory term if needed
  if (configuration.grab<bool>("with_memory"))
  {
    // build -lambda*a
    Vec_unique disp = create_unique_vec();
    perr = VecDuplicate(m_p_map->displacements(), disp.get());
    CHKERRABORT(m_comm, perr);
    perr = VecCopy(m_p_map->displacements(), *disp);
    CHKERRABORT(m_comm, perr);
    perr = VecScale(*disp, -total_mult);
    CHKERRABORT(m_comm, perr);

    // calculate vecdot(lapl_2, -lambda*a) and add to rhs in one operation
    perr = MatMultAdd(*m_p_map->laplacian(), *disp, *m_workspace->m_rhs, *m_workspace->m_rhs);
    CHKERRABORT(m_comm, perr);
  }

  // Force free tmat as no longer needed
  m_workspace->m_tmat = create_unique_mat();

  // solve for delta a
  KSP_unique m_ksp = create_unique_ksp();
  perr = KSPCreate(m_comm, m_ksp.get());
  CHKERRABORT(m_comm, perr);
  perr = KSPSetOperators(*m_ksp, *normmat, *normmat);
  CHKERRABORT(m_comm, perr);
  perr = KSPSetUp(*m_ksp);
  CHKERRABORT(m_comm, perr);
  perr = KSPSetFromOptions(*m_ksp);
  CHKERRABORT(m_comm, perr);
  perr = KSPSetUp(*m_ksp);
  CHKERRABORT(m_comm, perr);
  perr = KSPSolve(*m_ksp, *m_workspace->m_rhs, *m_workspace->m_delta);
  CHKERRABORT(m_comm, perr);
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
  auto start_iter = currspc.begin();
  auto end_iter = std::next(start_iter, m_fixed.ndim());
  while (all_true_varlen(start_iter, end_iter, imshape.begin(), imshape.end(),
      [](floating x, integer y) -> bool { return (y / x) > 2.0; }))
  {
    std::transform(start_iter, end_iter, start_iter, [](floating a) -> floating { return a * 2; });
    m_v_nodespacings.push_back(currspc);
  }
}
// iternum may be unused depending on debug level
void Elastic::calculate_tmat(integer iternum __attribute__((unused)))
{
  // Calculate average intensity 0.5(f+m)
  // Constant offset needed later, does not affect gradients
  PetscErrorCode perr = VecSet(*m_workspace->m_globaltmps[m_fixed.ndim()], -1.0);
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  // NB Z = aX + bY + cZ has call signature VecAXPBYPCZ(Z, a, b, c, X, Y) because reasons....
  perr = VecAXPBYPCZ(*m_workspace->m_globaltmps[m_fixed.ndim()], 0.5, 0.5, 1,
      *m_fixed.global_vec(), *m_p_registered->global_vec());
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter this to local for later
  perr = DMGlobalToLocalBegin(*m_fixed.dmda(), *m_workspace->m_globaltmps[m_fixed.ndim()],
      INSERT_VALUES, *m_workspace->m_localtmp);
  CHKERRABORT(m_comm, perr);
  perr = DMGlobalToLocalEnd(*m_fixed.dmda(), *m_workspace->m_globaltmps[m_fixed.ndim()],
      INSERT_VALUES, *m_workspace->m_localtmp);
  CHKERRABORT(m_comm, perr);

  // find average gradients
  for (uinteger idim = 0; idim < m_fixed.ndim(); idim++)
  {
    // likely change this to avoid mallocs/frees
    fd::gradient_existing(
        *(m_fixed.dmda()), *m_workspace->m_localtmp, *m_workspace->m_globaltmps[idim], idim);
  }

  // Negate average intensity to get 1 - 0.5(f+m) as needed by algorithm
  perr = VecScale(*m_workspace->m_globaltmps[m_fixed.ndim()], -1.0);
  CHKERRABORT(PETSC_COMM_WORLD, perr);

  // scatter grads into stacked vector
  m_workspace->scatter_grads_to_stacked();

  // 3. copy basis into p_tmat
  m_workspace->m_tmat = create_unique_mat();
  perr = MatDuplicate(*m_p_map->basis(), MAT_COPY_VALUES, m_workspace->m_tmat.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*m_workspace->m_tmat, std::string("Mat_tmat_") + std::to_string(iternum));

  // 4. left diagonal multiply p_tmat with stacked vector
  perr = MatDiagonalScale(*m_workspace->m_tmat, *m_workspace->m_stacktmp, nullptr);
  CHKERRABORT(m_comm, perr);
}

void Elastic::save_debug_frame(integer outer_count, integer inner_count)
{
  std::ostringstream outname;
  std::string file_str(configuration.grab<std::string>("intermediate_template"));
  bf::path registered_path(configuration.grab<std::string>("registered"));

  boost::format pad2("%02d");
  replace_token(file_str, ConfigurationBase::k_outer_token, (pad2 % outer_count).str());

  boost::format pad3("%03d");
  replace_token(file_str, ConfigurationBase::k_inner_token, (pad3 % inner_count).str());

  replace_token(file_str, ConfigurationBase::k_stem_token, 
                registered_path.filename().stem().string());

  replace_token(file_str, ConfigurationBase::k_extension_token,
                registered_path.extension().string());

  bf::path output_path(registered_path.parent_path());

  bf::path intermediates_path(configuration.grab<std::string>("intermediate_directory"));
  if(intermediates_path.is_absolute())
  {
    if (!bf::exists(intermediates_path))
    {
      std::ostringstream errss;
      errss << "Intermediate frame output path " << intermediates_path << " does not exist.";
      throw std::runtime_error(errss.str());
    }
    output_path = intermediates_path;
  }
  else
  {
    // Don't create directories we don't say we will, error instead
    // (Pre-flight checks mean this should never throw)
    if (!output_path.empty() && !bf::exists(output_path))
    {
      std::ostringstream errss;
      errss << "Output path " << output_path << " does not exist.";
      throw std::runtime_error(errss.str());
    }

    // Append intermediates dir
    output_path /= configuration.grab<std::string>("intermediate_directory");

    // Create intermediates dir as needed
    bf::create_directories(output_path);
  }

  // Finally append filename and any group extension
  output_path /= file_str;
  std::string output_path_str(output_path.string());
  output_path_str.append(":");
  output_path_str.append(configuration.grab<std::string>("registered_h5_path"));

  BaseWriter_unique wtr = BaseWriter::get_writer_for_filename(output_path_str, m_comm);
  wtr->write_image(*m_p_registered);
}

floating Elastic::approximate_optimum_lambda(Mat& mat_a, Mat& mat_b, floating lambda_mult,
    floating initial_guess, floating search_width, uinteger max_iter, floating lambda_min)
{
  floating x_lo =
      initial_guess - search_width > lambda_min ? initial_guess - search_width : lambda_min;
  floating x_mid = initial_guess;
  floating x_hi = initial_guess + search_width;

  floating y_lo(0), y_mid(0), y_hi(0);
  bool recalc_first = true;

  for (uinteger iter = 0; iter < max_iter; iter++)
  {
    Mat_unique mat_c = create_unique_mat();
    PetscErrorCode perr = MatDuplicate(mat_a, MAT_COPY_VALUES, mat_c.get());

    // calculate at lower value
    if (recalc_first)
    {
      perr = MatAXPY(*mat_c, lambda_mult * x_lo, mat_b, DIFFERENT_NONZERO_PATTERN);
      CHKERRABORT(m_comm, perr);
      y_lo = get_condnum_by_poweriter(*mat_c, 0.01, 100);
    }
    recalc_first = true;

    // Calculate at initial guess
    perr = MatAXPY(*mat_c, lambda_mult * (x_mid - x_lo), mat_b, DIFFERENT_NONZERO_PATTERN);
    CHKERRABORT(m_comm, perr);
    y_mid = get_condnum_by_poweriter(*mat_c, 0.01, 100);

    // calculate at higher value
    perr = MatAXPY(*mat_c, lambda_mult * (x_hi - x_mid), mat_b, DIFFERENT_NONZERO_PATTERN);
    CHKERRABORT(m_comm, perr);
    y_hi = get_condnum_by_poweriter(*mat_c, 0.01, 100);

    PetscPrintf(
        m_comm, " lo: %g(%g) mid: %g(%g) hi: %g(%g)\n", x_lo, y_lo, x_mid, y_mid, x_hi, y_hi);

    // check initial lower than either side
    if (y_lo >= y_mid)
    {
      // two options, either close to a minimum or need to search to higher x
      if (y_hi < y_mid)
      {
        // All points are negative gradient: keep searching to higher values of x:
        x_mid *= 2;
        x_hi *= 2;
        recalc_first = false;
        continue;
      }
      // Otherwise we have are close to a local minimum and are done
      break;
    }
    else
    {
      if (x_lo == lambda_min)
      {
        break;
      }
      // All negative gradient so keep searching left
      x_lo = x_lo >= 2.0 * lambda_min ? x_lo / 2 : 1.0 * lambda_min;
      x_mid = x_mid >= 2.0 * x_lo ? x_mid / 2 : x_lo + 1.0;
      x_hi = x_hi >= 2.0 * x_mid ? x_hi / 2 : x_mid + 1.0;
      continue;
    }
  }
  // fit quadratic and return minimum
  floating a, b, c;
  quadratic_from_points(x_lo, x_mid, x_hi, y_lo, y_mid, y_hi, a, b, c);

  floating x, y;
  quadratic_vertex(a, b, c, x, y);

  if (x < lambda_min)
  {
    x = lambda_min;
  }

  return x;
}
