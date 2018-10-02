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

  Map(const Image& mask, const floatvector node_spacing);
  Map(const Map& map, const floatvector node_spacing);

//  ~Map();

  Mat* basis() { return m_basis.get();}

  std::unique_ptr<Map> interpolate(floatvector new_spacing);

  //private:

  MPI_Comm m_comm;
  integer m_ndim;
  floatvector m_v_node_spacing;
  floatvector m_v_offsets;
  intvector m_v_image_shape;
  intvector m_v_map_shape;
  floatvector2d m_vv_node_locs;
  Mat_unique m_basis;
  Vec_unique m_displacements;

  void alloc_displacements();
  void calculate_node_locs();
  void calculate_basis();

};

#endif
