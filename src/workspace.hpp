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

#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include "elastic.hpp"
#include "image.hpp"
#include "types.hpp"

class WorkSpace {
public:
  WorkSpace(const Image& image, const Map& map);

  void reallocate_ephemeral_workspace(const Map& map);
  void scatter_stacked_to_grads();
  void scatter_grads_to_stacked();
  void duplicate_single_grad_to_stacked(size_t idx);

  void scatter_stacked_to_grads_noreorder();

  friend Elastic;
  friend Map;

  //  protected:

  void allocate_persistent_workspace();
  void create_reordering_scatterers();
  void create_nonreordering_scatterers();

  MPI_Comm m_comm;
  DM_shared m_dmda;
  integer m_size;
  std::vector<Vec_unique> m_globaltmps;
  std::vector<IS_unique> m_iss;
  std::vector<VecScatter_unique> m_r_scatterers;
  std::vector<VecScatter_unique> m_nr_scatterers;
  Vec_unique m_stacktmp, m_localtmp;
  Vec_unique m_delta, m_rhs;
  Mat_unique m_tmat;

  integer ephemeral_count;
};

#endif
