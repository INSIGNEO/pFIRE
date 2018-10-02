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

    void autoregister();


//  protected:

    integer m_max_iter = 50;

    // Straightforward initialize-by-copy
    MPI_Comm m_comm;
    integer m_imgdims;
    integer m_mapdims;
    integer m_size;
    integer m_iternum;
    const Image& m_fixed;
    const Image& m_moved;

    // Other class data, default initialize then populate in c'tor
    floatvector m_v_final_nodespacing;
    floatvector2d m_v_nodespacings;
    std::shared_ptr<Image> m_p_registered;
    std::vector<Vec_unique> m_vp_grads;
    Vec_unique m_p_stacked_grads;
    Mat_unique m_p_tmat;
    std::unique_ptr<Map> m_p_map;
    std::vector<IS_unique> m_vp_iss;
    std::vector<VecScatter_unique> m_vp_scatterers;

    void innerloop();
    void innerstep();

    void calculate_node_spacings();
    void allocate_scratch_storage();
    void create_scatterers();
    
    void scatter_grads_to_stacked();
    void create_t_matrix();
};

#endif
