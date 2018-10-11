#include "map.hpp"

Map::Map(const Image& mask, const floatvector node_spacing)
  : m_mask(mask), m_comm(mask.comm()), m_ndim(mask.ndim()), m_v_node_spacing(node_spacing),
    m_v_image_shape(mask.shape()), m_v_map_shape(intvector()), m_vv_node_locs(floatvector2d()),
    m_displacements(create_unique_vec())
{
  calculate_node_locs();
  calculate_basis();
  // TODO: mask basis
  // initialize displacement storage
  alloc_displacements();
  calculate_laplacian();
}

Map::Map(const Map& map, const floatvector new_spacing)
  : m_comm(map.m_comm), m_mask(map.m_mask), m_ndim(map.m_ndim), m_v_node_spacing(new_spacing),
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

void Map::update(const Vec &delta_vec)
{
  PetscErrorCode perr = VecAXPY(*m_displacements, -1, delta_vec);CHKERRABORT(m_comm, perr);
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

std::unique_ptr<Image> Map::warp(const Image& image, WorkSpace& wksp)
{
  // TODO: Check image is compatible
  
  // interpolate map to image nodes with basis
  PetscErrorCode perr = MatMult(*m_basis, *m_displacements,
                                *wksp.m_stacktmp);CHKERRABORT(m_comm, perr);
  wksp.scatter_stacked_to_grads();
  
  // build warp matrix
  std::vector<Vec*> tmps(0);
  for (auto const& vptr: wksp.m_globaltmps){ tmps.push_back(vptr.get());}
  Mat_unique warp = build_warp_matrix(m_comm, m_v_image_shape, tmps); 

  // apply matrix to get new image data
  std::unique_ptr<Image> new_img = image.duplicate();
  perr = MatMult(*warp, *image.global_vec(), *new_img->global_vec());CHKERRABORT(m_comm, perr);
  
  return new_img;
}

void Map::calculate_basis()
{
  // Get the full Nd basis
  m_basis = build_basis_matrix(m_comm, m_v_map_shape, m_v_image_shape,
                               m_v_node_spacing, m_v_offsets, m_ndim+1);

  // Now grab a 1d basis as a submatrix. Note can't do this the other way round because Petsc won't
  // allow reuse of rows/cols in MatCreateSubMatrix
  // Work out what rows, columns the rank needs to own for compatability with image and
  // displacement vectors
/*  integer rowstart, rowsize, colstart, colsize;
  PetscErrorCode perr = VecGetOwnershipRange(*m_mask.global_vec(), &rowstart, &rowsize);CHKERRABORT(m_comm, perr);
  rowsize -= rowstart;
  throw std::runtime_error("need to get just one map length's worth of columns");
  perr = MatGetOwnershipRangeColumn(*m_basis, &colstart, &colsize);CHKERRABORT(m_comm, perr);
  colsize -= colstart;
  
  // Express these ranges as index sets 
  IS_unique rows = create_unique_is();
  perr = ISCreateStride(m_comm, rowsize, rowstart, 1, rows.get());CHKERRABORT(m_comm, perr);
  IS_unique cols = create_unique_is();
  perr = ISCreateStride(m_comm, colsize, colstart, 1, cols.get());CHKERRABORT(m_comm, perr);

  m_basis_1d = create_unique_mat();
  perr = MatCreateSubMatrix(*m_basis, *rows, *cols, MAT_INITIAL_MATRIX, m_basis_1d.get());*/
}

void Map::calculate_laplacian()
{
  integer startrow, endrow;
  PetscErrorCode perr = VecGetOwnershipRange(*m_displacements, &startrow, &endrow);CHKERRABORT(m_comm, perr);
  m_lapl = build_laplacian_matrix(m_comm, m_v_map_shape, startrow, endrow, m_ndim+1);
}
