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

#include "baseconfiguration.hpp"
#include "types.hpp"

class ElasticRegistration {
public:
  ElasticRegistration(std::shared_ptr<Image> fixed, std::shared_ptr<Image> moved, std::shared_ptr<Mask> mask,
                      const intcoord& target_spacing, const ConfigurationBase& config);

  void set_initial_guess(const MapBase& initial_guess);

  void autoregister(integer max_iterations = 50);

  std::shared_ptr<Image> registered() const { return _registered; }

  const MapBase& result_map() { return *_current_map;}

protected:
  void registration_inner_loop(integer max_iterations);
  Vec_unique solver_step(floating internum, bool recalculate_lambda);
  bool is_converged(const Vec& map_delta, floating threshold = 0.1);

  void save_debug_frame(integer iteration_num);
  floating approximate_optimum_lambda(const Mat& mat_a, const Mat& mat_b, floating lambda_mult,
                                      floating initial_guess, floating search_width, uinteger max_iter,
                                      floating lambda_min);

private:
  const MPI_Comm _comm;
  std::shared_ptr<Image> _fixed;
  std::shared_ptr<Image> _moved;
  std::shared_ptr<Mask> _mask;
  intcoord _target_spacing;
  const ConfigurationBase& _configuration;

  std::unique_ptr<MapManager> _mapmanager;
  std::shared_ptr<MapBase> _current_map;
  std::shared_ptr<Image> _registered;

  const size_t _mi_window = 5;
  std::deque<floating> _previous_mis;

  floating _lambda = 10;
  floating _initial_search_width = 10;
  floating _lambda_search_maxiter = 5;
  floating _lambda_min = 1;
};
