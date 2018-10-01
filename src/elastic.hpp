#ifndef ELASTIC_HPP
#define ELASTIC_HPP

#include<algorithm>

#include "types.hpp"
#include "image.hpp"
#include "map.hpp"

class Elastic
{
  public:
    Elastic(const Image& fixed, const Image& moved, Map& map);

    void do_scatter_to_stacked();

  protected:
    MPI_Comm m_comm;
    integer m_ndim;
    integer m_size;
    std::vector<Vec_unique> m_vp_grads;
    Vec_unique m_p_stacked_grads;
    std::vector<IS_unique> m_vp_iss;
    std::vector<VecScatter_unique> m_vp_scatterers;

    void create_scatterers();
    
    void create_t_matrix()
};

#endif
