#include "map.hpp"

Map::Map(const Image& mask, const floatvector node_spacing)
  : m_comm(mask.comm()), m_ndim(mask.ndim()), m_v_node_spacing(node_spacing),
    m_v_image_shape(mask.shape()), m_v_map_shape(intvector()), m_vv_node_locs(floatvector2d()),
    m_displacements(create_unique_vec())
{
  std::cout << "Begin map ctr" << std::endl;
  calculate_node_locs();
  std::cout << "Node locs calculated" << std::endl;
  calculate_basis();
  std::cout << "Basis calculated" << std::endl;
  // TODO: mask basis
  // initialize displacement storage
  alloc_displacements();
  std::cout << "Memory allocated" << std::endl;
  calculate_laplacian();
  std::cout << "Laplacian calculated" << std::endl;
  std::cout << "End map ctr" << std::endl;
}

Map::Map(const Map& map, const floatvector new_spacing)
  : m_comm(map.m_comm), m_ndim(map.m_ndim), m_v_node_spacing(new_spacing),
    m_v_image_shape(map.m_v_image_shape), m_v_map_shape(intvector()), 
    m_vv_node_locs(floatvector2d()), m_displacements(create_unique_vec())
{
  calculate_node_locs();
  calculate_basis();
  // TODO: mask basis
  // initialize displacement storage
  alloc_displacements();
  calculate_laplacian();
}

std::unique_ptr<Map> Map::interpolate(floatvector new_spacing)
{
  std::unique_ptr<Map> new_map(new Map(*this, new_spacing));

  floatvector scalings(m_ndim, 0.0);
  floatvector offsets(m_ndim, 0.0);

  n_ary_transform(std::divides<>(), scalings.begin(), this->m_v_node_spacing.begin(),
                  this->m_v_node_spacing.end(), new_map->m_v_node_spacing.begin());
  n_ary_transform(std::minus<>(), offsets.begin(), new_map->m_v_offsets.begin(),
                  new_map->m_v_offsets.end(), this->m_v_offsets.begin());

  Mat_unique interp = build_basis_matrix(m_comm, m_v_map_shape, new_map->m_v_map_shape,
                                         scalings, offsets, m_ndim+1);

  PetscErrorCode perr = MatMult(*interp, *m_displacements, *new_map->m_displacements);CHKERRABORT(m_comm, perr);

  return new_map;
}


void Map::alloc_displacements()
{
  MatCreateVecs(*m_basis, m_displacements.get(), nullptr);  
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
    floatvector nodes(num_nod, 0.0);
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

void Map::calculate_laplacian()
{
  std::cout << "Begin lapl calc" << std::endl;
  integer startrow, endrow;
  PetscErrorCode perr = VecGetOwnershipRange(*m_displacements, &startrow, &endrow);CHKERRABORT(m_comm, perr);
  m_lapl = build_laplacian_matrix(m_comm, m_v_map_shape, startrow, endrow, m_ndim+1);
  std::cout << "End lapl calc" << std::endl;
}
