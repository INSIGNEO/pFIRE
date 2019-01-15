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

  friend Elastic;
  friend Map;

  //  protected:

  void allocate_persistent_workspace();
  void create_scatterers();

  MPI_Comm m_comm;
  DM_shared m_dmda;
  integer m_size;
  std::vector<Vec_unique> m_globaltmps;
  std::vector<IS_unique> m_iss;
  std::vector<VecScatter_unique> m_scatterers;
  Vec_unique m_stacktmp, m_localtmp;
  Vec_unique m_delta, m_rhs;
  Mat_unique m_tmat;

  integer ephemeral_count;
};

#endif
