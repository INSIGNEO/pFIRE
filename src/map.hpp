#ifndef MAP_HPP
#define MAP_HPP

#include<exception>
#include<numeric>
#include<iostream>

#include "types.hpp"

class Map{

  public:

  Map(const Image& mask, const floatvector& node_spacing);
  Map(const Map& map, const floatvector& node_spacing);

//  ~Map();

  Mat* basis() const{ return m_basis.get();}
  Mat* laplacian() const{ return m_lapl.get();}
  const floatvector2d node_locs() const{ return m_vv_node_locs;}

  const MPI_Comm& comm() const{ return m_comm;}
  uinteger ndim() const{ return m_ndim;}
  const intvector& shape() const{ return m_v_map_shape;}
  integer size() const
  { 
    return std::accumulate(m_v_map_shape.cbegin(), m_v_map_shape.cend(), 1, std::multiplies<>());
  }

  std::pair<integer, integer> get_displacement_ownershiprange() const;

  const floating* get_raw_data_ro() const;
  void release_raw_data_ro(const floating*& ptr) const;

  void update(const Vec &delta_vec);
  std::unique_ptr<Map> interpolate(const floatvector& new_spacing);

  std::unique_ptr<Image> warp(const Image& image, WorkSpace& wksp);

  //private:

  MPI_Comm m_comm;
  const Image& m_mask;
  uinteger m_ndim;
  floatvector m_v_node_spacing;
  floatvector m_v_offsets;
  intvector m_v_image_shape;
  intvector m_v_map_shape;
  floatvector2d m_vv_node_locs;
  Mat_unique m_basis;
  Mat_unique m_lapl;
  Vec_unique m_displacements;

  void alloc_displacements();
  void calculate_node_locs();
  void calculate_basis();
  void calculate_laplacian();
  void calculate_warp_matrix();

};

#endif
