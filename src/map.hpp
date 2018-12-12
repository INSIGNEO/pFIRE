#ifndef MAP_HPP
#define MAP_HPP

#include <exception>
#include <iostream>
#include <numeric>
#include <utility>

#include "types.hpp"

class Map {
public:
  Map(const Image& mask, const floatvector& node_spacing);
  Map(const Map& map, const floatvector& node_spacing);

  //  ~Map();

  Mat* basis() const
  {
    return m_basis.get();
  }
  Mat* laplacian() const
  {
    return m_lapl.get();
  }
  const floatvector2d node_locs() const
  {
    return m_vv_node_locs;
  }
  const MPI_Comm& comm() const
  {
    return m_comm;
  }
  uinteger ndim() const
  {
    return m_ndim;
  }
  const intvector& shape() const
  {
    return map_shape;
  }
  const intvector& image_shape() const
  {
    return m_v_image_shape;
  }
  const floatvector& spacing() const
  {
    return m_v_node_spacing;
  }
  integer size() const
  {
    return std::accumulate(map_shape.cbegin(), map_shape.cend(), 1, std::multiplies<>());
  }

  floatvector low_corner() const;
  std::pair<integer, integer> get_displacement_ownershiprange() const;

  const floating* get_raw_data_ro() const;
  void release_raw_data_ro(const floating*& ptr) const;

  void update(const Vec& delta_vec);
  std::unique_ptr<Map> interpolate(const floatvector& new_spacing);

  std::unique_ptr<Image> warp(const Image& image, WorkSpace& wksp);

  std::pair<intvector, intvector> get_dmda_local_extents() const;
  Vec_unique get_dim_data_dmda_blocked(integer dim) const;

  static intvector
  calculate_map_shape(intvector const& image_shape, floatvector const& nodespacing);

  // private:

  MPI_Comm m_comm;
  const Image& m_mask;
  uinteger m_ndim;
  floatvector m_v_node_spacing;
  floatvector m_v_offsets;
  intvector m_v_image_shape;
  intvector map_shape;
  floatvector2d m_vv_node_locs;
  Mat_unique m_basis;
  Mat_unique m_lapl;
  Vec_unique m_displacements;
  mutable DM_unique map_dmda;

  void apply_mask_to_basis();
  void alloc_displacements();
  void initialize_dmda() const;
  void calculate_node_locs();
  void calculate_basis();
  void calculate_laplacian();
  void calculate_warp_matrix();
};

#endif
