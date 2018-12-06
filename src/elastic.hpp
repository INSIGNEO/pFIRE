#ifndef ELASTIC_HPP
#define ELASTIC_HPP

#include <algorithm>
#include <iostream>

#include <petscdmda.h>
#include <petscmat.h>

#include "image.hpp"
#include "map.hpp"
#include "types.hpp"
#include "workspace.hpp"

class Elastic {
public:
  Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing);

  void autoregister();

  std::shared_ptr<Image> registered() const
  {
    return m_p_registered;
  }

  //  protected:

  integer m_max_iter = 50;
  floating m_convergence_thres = 0.1;

  // Straightforward initialize-by-copy
  MPI_Comm m_comm;
  integer m_imgdims;
  integer m_mapdims;
  integer m_size;
  integer m_iternum;
  const Image& m_fixed;
  const Image& m_moved;

  // Other class data, default initialize then populate in c'tor
  floatvector2d m_v_nodespacings;
  floatvector m_v_final_nodespacing;
  std::shared_ptr<Image> m_p_registered;
  std::unique_ptr<Map> m_p_map;
  std::shared_ptr<WorkSpace> m_workspace;
  Mat_unique normmat;

  void save_debug_frame(integer ocount, integer icount);
  void innerloop(integer outer_count);
  void innerstep(floating lambda);

  void block_precondition();
  void calculate_node_spacings();
  void calculate_tmat();
};

#endif
