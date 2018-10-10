#ifndef ELASTIC_HPP
#define ELASTIC_HPP

#include<algorithm>
#include<iostream>

#include<petscmat.h>
#include<petscdmda.h>

#include "types.hpp"
#include "image.hpp"
#include "workspace.hpp"
#include "map.hpp"

class Elastic
{
  public:
    Elastic(const Image& fixed, const Image& moved, const floatvector nodespacing);

    void autoregister();


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
    floatvector m_v_final_nodespacing;
    floatvector2d m_v_nodespacings;
    std::shared_ptr<Image> m_p_registered;
    std::unique_ptr<Map> m_p_map;
    std::shared_ptr<WorkSpace> m_workspace;

    Mat_unique m_p_normmat;

    void innerloop();
    void innerstep(floating lambda);

    void calculate_node_spacings();
    void calculate_tmat();
};

#endif
