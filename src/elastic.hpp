#ifndef ELASTIC_HPP
#define ELASTIC_HPP

#include<algorithm>

#include<petscmat.h>

#include "types.hpp"
#include "image.hpp"
#include "map.hpp"

class Elastic
{
  public:
    Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing);

    void do_scatter_to_stacked();

  protected:
    MPI_Comm m_comm;
    integer m_imgdims;
    integer m_mapdims;
    integer m_size;
    const Image& m_fixed;
    const Image& m_moved;
    std::vector<Vec_unique> m_vp_grads;
    Vec_unique m_p_stacked_grads;
    std::unique_ptr<Map> m_p_map;
    std::vector<IS_unique> m_vp_iss;
    std::vector<VecScatter_unique> m_vp_scatterers;


    void create_scratch();
    void create_scatterers();
    
    void create_t_matrix()
};

#endif
