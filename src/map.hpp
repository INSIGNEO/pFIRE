#ifndef MAP_HPP
#define MAP_HPP

#include<exception>
#include<numeric>
#include<iostream>

#include "types.hpp"
#include "indexing.hpp"
#include "image.hpp"
#include "basis.hpp"
#include "laplacian.hpp"
#include "workspace.hpp"

class Map{

  public:

  Map(const Image& mask, const floatvector& node_spacing);
  Map(const Map& map, const floatvector& node_spacing);

//  ~Map();

  Mat* basis() const{ return m_basis.get();}
  Mat* laplacian() const{ return m_lapl.get();}

  integer size() const
  { return std::accumulate(m_v_map_shape.cbegin(), m_v_map_shape.cend(), 1, std::multiplies<>());
  }

  void update(const Vec &delta_vec);
  std::unique_ptr<Map> interpolate(const floatvector& new_spacing);

  std::unique_ptr<Image> warp(const Image& image, WorkSpace& wksp);

  //private:

  MPI_Comm m_comm;
  const Image& m_mask;
  integer m_ndim;
  floatvector m_v_node_spacing;
  floatvector m_v_offsets;
  intvector m_v_image_shape;
  intvector m_v_map_shape;
  floatvector2d m_vv_node_locs;
  Mat_unique m_basis;
//  Mat_unique m_basis_1d;
  Mat_unique m_lapl;
  Vec_unique m_displacements;

  void alloc_displacements();
  void calculate_node_locs();
  void calculate_basis();
  void calculate_laplacian();
  void calculate_warp_matrix();

};

#endif
