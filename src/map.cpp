#include "map.hpp"

Map::Map(const intvector node_spacing, const Image& mask)
  : m_comm(mask.comm()), m_ndim(mask.ndim()), m_v_node_spacing(node_spacing),
    m_v_image_shape(mask.shape()), m_v_map_shape(intvector()), m_vv_node_locs(floatvector2d()),
    m_dmda(mask.dmda())
{
  calculate_node_locs();
  calculate_basis();
  // TODO: mask basis
  // initialize displacement storage
}

void Map::calculate_node_locs()
{
  // calculate self size and offset
  for(integer idim=0; idim < m_ndim; idim++)
  {
    floating centre = m_v_image_shape[idim]/2. - 0.5;
    // want always to have odd number of nodes so find num nodes for each half,
    // multiply by two and subtract one to get total nodes
    // N.B number of nodes = 1 + number of spaces, so:
    integer num_spc = (integer)std::ceil(centre/m_v_node_spacing[idim]);
    integer num_nod = num_spc*2 + 1;
    floating lo = centre - num_spc*m_v_node_spacing[idim];
    floatvector nodes(num_nod);
    std::generate(nodes.begin(), nodes.end(),
                  [x = lo, y = m_v_node_spacing[idim]] () mutable { x += y; return x;});
    m_vv_node_locs.push_back(nodes);
    m_v_map_shape.push_back(num_nod);
    m_v_offsets.push_back(lo);
  }
}


void Map::calculate_basis()
{
  m_basis = build_basis_matrix(m_comm, m_v_map_shape, m_v_image_shape,
                               m_v_node_spacing, m_v_offsets, m_ndim+1);
}

