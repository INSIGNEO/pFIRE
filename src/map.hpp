#ifndef MAP_HPP
#define MAP_HPP

#include<numeric>
#include<iostream>

#include "types.hpp"
#include "indexing.hpp"
#include "image.hpp"
#include "basis.hpp"

class Image;

class Map{

  public:

  Map(const intvector node_spacing, const Image& image);

//  ~Map();

  Mat* basis() { return m_basis.get();}

  Map interpolate(intvector new_spacing);

  private:

  MPI_Comm m_comm;
  integer m_ndim;
  intvector m_v_node_spacing;
  floatvector m_v_offsets;
  intvector m_v_image_shape;
  intvector m_v_map_shape;
  floatvector2d m_vv_node_locs;
  DM_shared m_dmda;
  Mat_unique m_basis;
  Vec_unique m_displacements;

  void calculate_node_locs();
  void calculate_basis();

};

#endif
