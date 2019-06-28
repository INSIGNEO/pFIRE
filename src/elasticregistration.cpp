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

#include "elasticregistration.hpp"

#include <boost/format.hpp>

#include "basewriter.hpp"
#include "exceptions.hpp"
#include "file_utils.hpp"
#include "image.hpp"
#include "mapmanager.hpp"
#include "mask.hpp"
#include "math_utils.hpp"
#include "tmatrix.hpp"

ElasticRegistration::ElasticRegistration(std::shared_ptr<Image> fixed, std::shared_ptr<Image> moved,
                                         std::shared_ptr<Mask> mask, const intcoord& target_spacing,
                                         const ConfigurationBase& config)
  : _comm(fixed->comm()), _fixed(std::move(fixed)), _moved(std::move(moved)), _mask(std::move(mask)),
    _target_spacing(target_spacing), _configuration(config),
    _mapmanager(std::make_unique<MapManager>(target_spacing, *_mask)),
    _registered(GridVariable::copy(*_moved))
{
}

void ElasticRegistration::set_initial_guess(const MapBase& initial_guess __attribute__((unused)))
{
  // Query manager for map closest in shape to initial_guess
  // Set Mapmanager index to that shape
  // Interpolate initial_guess to current map
  // Carry on from there...
  throw InternalError("Not yet implemented");
}

void ElasticRegistration::autoregister(integer max_iterations)
{
  for (_current_map = _mapmanager->get_current_map(); _current_map != nullptr;
       _current_map = _mapmanager->get_next_map())
  {
    registration_inner_loop(max_iterations);
  }
  _current_map = _mapmanager->get_current_map(); // Replace nullptr with final valid map
}

void ElasticRegistration::registration_inner_loop(integer max_iterations)
{
  integer iteration = 0;
  if (_configuration.grab<bool>("save_intermediate_frames"))
  {
    save_debug_frame(iteration);
  }
  for (iteration = 1; iteration <= max_iterations; iteration++)
  {
    this->_previous_mis.push_back(this->_fixed->mutual_information(*this->_moved));

    // Do solver step, recalculating lambda if it is the first iteration of the loop
    Vec_unique map_delta = solver_step(iteration, iteration == 1);

    if (this->is_converged(*map_delta))
    {
      break;
    }
  }
  if (iteration == max_iterations)
  {
    PetscPrintf(_comm, "Warning: Maximum iteration limit reached (%d)\n", max_iterations);
  }
}

Vec_unique ElasticRegistration::solver_step(floating iternum, bool recalculate_lambda)
{
  auto tmat2_and_tmatfm = build_tmat2_and_tmatfm(*this->_fixed, *_registered, *this->_current_map);
  Mat_unique tmat2 = std::move(tmat2_and_tmatfm.first);
  Vec_unique tmatfm = std::move(tmat2_and_tmatfm.second);

  floating lapl2_diag_sum = diagonal_sum(_current_map->laplacian2());
  floating tmat2_diag_sum = diagonal_sum(*tmat2);
  floating lapl_premult = tmat2_diag_sum / lapl2_diag_sum;
  PetscPrintf(_comm, "lapl_mult = %f, laplsum= %f, tmatsum = %f\n", lapl_premult, lapl2_diag_sum,
              tmat2_diag_sum);


  // Recalculate lambda if requested
  if (recalculate_lambda)
  {
    _lambda = approximate_optimum_lambda(*tmat2, _current_map->laplacian2(), lapl_premult, _lambda,
                                         _initial_search_width, _lambda_search_maxiter, _lambda_min);
    PetscPrintf(_comm, "Calculated lambda = %f\n", _lambda);
  }

  // Add lambda*L^tL to T^t T
  PetscErrorCode perr = MatAXPY(*tmat2, _lambda * lapl_premult, _current_map->laplacian2(),
                                DIFFERENT_NONZERO_PATTERN);
  CHKERRXX(perr);

  // Subtract lambda*L^tL*a from T^t(f-m)
  Vec_unique lambda_l2_a = create_unique_vec();
  perr = VecDuplicate(*tmatfm, lambda_l2_a.get());
  CHKERRXX(perr);
  perr = MatMult(_current_map->laplacian2(), _current_map->displacement_vector(), *lambda_l2_a);
  CHKERRXX(perr);
  //  perr = VecAXPY(*tmatfm, -_lambda * lapl_premult, *lambda_l2_a);
  // CHKERRXX(perr);

  // Co-opt now defunct lambda_l2_a as solution vector
  Vec_unique delta = std::move(lambda_l2_a);

  // Now do actual solve
  KSP_unique m_ksp = create_unique_ksp();
  perr = KSPCreate(_comm, m_ksp.get());
  CHKERRXX(perr);
  perr = KSPSetOperators(*m_ksp, *tmat2, *tmat2);
  CHKERRXX(perr);
  perr = KSPSetUp(*m_ksp);
  CHKERRXX(perr);
  perr = KSPSetFromOptions(*m_ksp);
  CHKERRXX(perr);
  perr = KSPSetUp(*m_ksp);
  CHKERRXX(perr);
  perr = KSPSolve(*m_ksp, *tmatfm, *delta);
  CHKERRXX(perr);

  // update map
  _current_map->update(*delta);

  // warp image
  _registered = _moved->warp(*_current_map);
  _registered->normalize();

  floating mutinf = _registered->mutual_information(*_fixed);
  _previous_mis.push_back(mutinf);
  PetscPrintf(_comm, "Mutual information: %f\n", mutinf);

  if (_configuration.grab<bool>("save_intermediate_frames"))
  {
    save_debug_frame(iternum);
  }

  return delta;
}

bool ElasticRegistration::is_converged(const Vec& map_delta, floating threshold)
{
  // Roll our own because we need to avoid intensity values
  std::array<floating, 2> avg_data;
  auto& avg_sum = avg_data[0];
  auto& avg_num = avg_data[1];

  floating localmax = PETSC_MIN_REAL;
  integer local_lo, local_hi;
  PetscErrorCode perr = VecGetOwnershipRange(map_delta, &local_lo, &local_hi);
  CHKERRXX(perr);
  avg_num = local_hi - local_lo;

  const floating* delta_data;
  perr = VecGetArrayRead(map_delta, &delta_data);
  CHKERRXX(perr);
  for (integer iidx = local_lo; iidx < local_hi; iidx++)
  {
    if (iidx % _current_map->ndof() == 0) // skip intensity values
    {
      continue;
    }
    floating absval = std::abs(delta_data[iidx]);
    localmax = absval > localmax ? absval : localmax;
    avg_sum += absval;
  }
  perr = VecRestoreArrayRead(map_delta, &delta_data);
  CHKERRXX(perr);

  MPI_Allreduce(MPI_IN_PLACE, &localmax, 1, MPIU_SCALAR, MPI_MAX, _comm);
  MPI_Allreduce(MPI_IN_PLACE, avg_data.data(), 2, MPIU_SCALAR, MPI_SUM, _comm);

  floating average = avg_sum / avg_num;

  PetscPrintf(_comm, "Displacement delta: %f / %f (avg / max)\n", average, localmax);

  bool disp_conv = average < threshold;

  // Now calculate gradient of MI over last n iterations
  // First check if we need to discard values
  while (this->_previous_mis.size() > this->_mi_window)
  {
    this->_previous_mis.pop_front();
  }

  // Now calculate simple linear regression
  intvector xpoints(this->_previous_mis.size());
  std::iota(xpoints.begin(), xpoints.end(), 1);
  floating num_points = static_cast<double>(xpoints.size());

  floating s_x = 0.5 * num_points * (num_points + 1);
  floating s_xx = (1.0 / 6.0) * (num_points * (num_points + 1) * (2 * num_points + 1));
  floating s_y = std::accumulate(this->_previous_mis.cbegin(), this->_previous_mis.cend(), 0.);
  floating s_xy = accumulate(xpoints.cbegin(), xpoints.cend(), this->_previous_mis.cbegin(), 0.,
                             [](floating x, integer y, floating z) -> floating { return x + y * z; });

  floating mi_grad = (num_points*s_xy - s_x*s_y)/(num_points*s_xx - s_x*s_x);

  PetscPrintf(_comm, "Mutual Information gradient: %f\n", mi_grad);

  bool mi_conv = false;
  if(this->_previous_mis.size() >= this->_mi_window)
  {
    mi_conv = mi_grad < 0;
  }

  return disp_conv || mi_conv;
}

void ElasticRegistration::save_debug_frame(integer iteration_num)
{
  integer outer_count = _mapmanager->step_number();
  std::ostringstream outname;
  std::string file_str(_configuration.grab<std::string>("intermediate_template"));
  bf::path registered_path(_configuration.grab<std::string>("registered"));

  boost::format pad2("%02d");
  replace_token(file_str, ConfigurationBase::k_outer_token, (pad2 % outer_count).str());

  boost::format pad3("%03d");
  replace_token(file_str, ConfigurationBase::k_inner_token, (pad3 % iteration_num).str());

  replace_token(file_str, ConfigurationBase::k_stem_token, registered_path.filename().stem().string());

  replace_token(file_str, ConfigurationBase::k_extension_token, registered_path.extension().string());

  bf::path output_path(registered_path.parent_path());

  bf::path intermediates_path(_configuration.grab<std::string>("intermediate_directory"));
  if (intermediates_path.is_absolute())
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
    output_path /= _configuration.grab<std::string>("intermediate_directory");

    // Create intermediates dir as needed
    bf::create_directories(output_path);
  }

  // Finally append filename and any group extension
  output_path /= file_str;
  std::string output_path_str(output_path.string());
  output_path_str.append(":");
  output_path_str.append(_configuration.grab<std::string>("registered_h5_path"));

  BaseWriter_unique wtr = BaseWriter::get_writer_for_filename(output_path_str, _comm);
  wtr->write_image(*_registered);
}

floating ElasticRegistration::approximate_optimum_lambda(const Mat& mat_a, const Mat& mat_b,
                                                         floating lambda_mult, floating initial_guess,
                                                         floating search_width, uinteger max_iter,
                                                         floating lambda_min)
{
  floating x_lo = initial_guess - search_width;
  x_lo = x_lo > lambda_min ? x_lo : lambda_min;
  floating x_mid = initial_guess;
  x_mid = x_mid > x_lo ? x_mid : x_lo + search_width/2;
  floating x_hi = initial_guess + search_width;
  x_hi = x_hi > x_mid ? x_hi : x_mid + search_width;

  floating y_lo(0), y_mid(0), y_hi(0);
  for (uinteger iter = 0; iter < max_iter; iter++)
  {
    Mat_unique mat_c = create_unique_mat();
    PetscErrorCode perr = MatDuplicate(mat_a, MAT_COPY_VALUES, mat_c.get());

    // calculate at lower value
    perr = MatAXPY(*mat_c, lambda_mult * x_lo, mat_b, DIFFERENT_NONZERO_PATTERN);
    CHKERRXX(perr);
    y_lo = get_condnum_by_poweriter(*mat_c, 0.01, 100);

    // Calculate at initial guess
    perr = MatAXPY(*mat_c, lambda_mult * x_mid, mat_b, DIFFERENT_NONZERO_PATTERN);
    CHKERRXX(perr);
    y_mid = get_condnum_by_poweriter(*mat_c, 0.01, 100);

    // calculate at higher value
    perr = MatAXPY(*mat_c, lambda_mult * x_hi, mat_b, DIFFERENT_NONZERO_PATTERN);
    CHKERRXX(perr);
    y_hi = get_condnum_by_poweriter(*mat_c, 0.01, 100);

    PetscPrintf(_comm, " lo: %g(%g) mid: %g(%g) hi: %g(%g)\n", x_lo, y_lo, x_mid, y_mid, x_hi, y_hi);

    // check initial lower than either side
    if (y_lo >= y_mid)
    {
      // two options, either close to a minimum or need to search to higher x
      if (y_hi < y_mid)
      {
        // All points are negative gradient: keep searching to higher values of x:
        x_mid *= 2;
        x_hi *= 2;
        continue;
      }
      // Otherwise we are close to a local minimum and are done
      break;
    }

    // If we already reached lambda_min we need to be done at that
    if (x_lo == lambda_min)
    {
      break;
    }

    // Otherwise it is a positive gradient so keep searching left
    x_lo = x_lo >= 2.0 * lambda_min ? x_lo / 2 : 1.0 * lambda_min;
    x_mid = x_mid >= 2.0 * x_lo ? x_mid / 2 : x_lo + 1.0;
    x_hi = x_hi >= 2.0 * x_mid ? x_hi / 2 : x_mid + 1.0;
  }
  // fit quadratic and return minimum
  floating a, b, c;
  quadratic_from_points(x_lo, x_mid, x_hi, y_lo, y_mid, y_hi, a, b, c);

  if (a < 0)
  {
    PetscPrintf(_comm, "Warning: convex function in smoothing parameter estimation.\n");
    return lambda_min;
  }

  floating x, y;
  quadratic_vertex(a, b, c, x, y);

  if (x < lambda_min)
  {
    return lambda_min;
  }

  return x;
}
