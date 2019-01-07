#include "map.hpp"

#include <petscvec.h>

#include "basis.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "laplacian.hpp"
#include "workspace.hpp"

#include "iterator_routines.hpp"

Map::Map(const Image& mask, const floatvector& node_spacing)
    : m_comm(mask.comm()), m_mask(mask), m_ndim(mask.ndim()), m_v_node_spacing(node_spacing),
      m_v_offsets(floatvector()), m_v_image_shape(mask.shape()), map_shape(intvector()),
      m_vv_node_locs(floatvector2d()), m_basis(create_unique_mat()), m_lapl(create_unique_mat()),
      m_displacements(create_unique_vec()), map_dmda(create_unique_dm())
{
  calculate_node_locs();
  calculate_basis();
  // TODO: mask basis
  // initialize displacement storage
  alloc_displacements();
  calculate_laplacian();
}

void Map::update(const Vec& delta_vec)
{
  PetscErrorCode perr = VecAXPY(*m_displacements, 1, delta_vec);
  CHKERRABORT(m_comm, perr);
}

std::unique_ptr<Map> Map::interpolate(const floatvector& new_spacing)
{
  std::unique_ptr<Map> new_map(new Map(this->m_mask, new_spacing));

  floatvector scalings(m_ndim, 0.0);
  floatvector offsets(m_ndim, 0.0);

  n_ary_transform(
      std::divides<>(), scalings.begin(), new_map->m_v_node_spacing.begin(),
      new_map->m_v_node_spacing.end(), this->m_v_node_spacing.begin());
  n_ary_transform(
      [](floating x, floating y, floating a) -> floating { return (x - y) / a; }, offsets.begin(),
      new_map->m_v_offsets.begin(), new_map->m_v_offsets.end(), this->m_v_offsets.begin(),
      this->m_v_node_spacing.begin());

  Mat_unique interp = build_basis_matrix(
      m_comm, map_shape, new_map->map_shape, scalings, offsets, m_ndim, m_ndim + 1);

  PetscErrorCode perr = MatMult(*interp, *m_displacements, *new_map->m_displacements);
  CHKERRABORT(m_comm, perr);

  return new_map;
}

void Map::initialize_dmda() const
{
  // Only initialize if needed
  if (*map_dmda != nullptr)
  {
    return;
  }

  integer dof_per_node = 1;
  integer stencil_width = 1;
  PetscErrorCode perr;

  // Make sure things get gracefully cleaned up
  perr = DMDACreate3d(
      m_comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, // BCs
      DMDA_STENCIL_STAR,                                            // stencil shape
      map_shape[0], map_shape[1], map_shape[2],                     // global mesh shape
      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,                     // ranks per dim
      dof_per_node, stencil_width,                                  // dof per node, stencil size
      nullptr, nullptr, nullptr, // partition sizes nullptr -> petsc chooses
      map_dmda.get());
  CHKERRABORT(m_comm, perr);

  perr = DMSetUp(*(map_dmda));
  CHKERRABORT(m_comm, perr);
}

void Map::alloc_displacements()
{
  MatCreateVecs(*m_basis, m_displacements.get(), nullptr);
}

intvector Map::calculate_map_shape(intvector const& image_shape, floatvector const& nodespacing)
{
  if (image_shape.size() != nodespacing.size())
  {
    throw std::runtime_error("Image and nodespacing dimensions must match");
  }

  // want always to have odd number of nodes so find num nodes for each half,
  // multiply by two and subtract one to get total nodes
  // N.B number of nodes = 1 + number of spaces, so:
  intvector nnod(image_shape.size());
  std::transform(
      image_shape.cbegin(), image_shape.cend(), nodespacing.cbegin(), nnod.begin(),
      [](integer is, floating ns) -> integer {
        return static_cast<integer>(1 + 2 * std::ceil((is / 2. - 0.5) / ns));
      });
  return nnod;
}

void Map::calculate_node_locs()
{
  map_shape = calculate_map_shape(m_v_image_shape, m_v_node_spacing);
  // calculate self size and offset
  for (uinteger idim = 0; idim < m_ndim; idim++)
  {
    floating lo = 0.5 * (m_v_image_shape[idim] - (map_shape[idim] - 1) * m_v_node_spacing[idim]);
    floatvector nodes(map_shape[idim], 0.0);
    std::generate(nodes.begin(), nodes.end(), [x = lo, y = m_v_node_spacing[idim]]() mutable {
      x += y;
      return x-y;
    });
    m_vv_node_locs.push_back(nodes);
    m_v_offsets.push_back(lo);
  }
}

std::unique_ptr<Image> Map::warp(const Image& image, WorkSpace& wksp)
{
  // TODO: Check image is compatible

  // interpolate map to image nodes with basis
  PetscErrorCode perr = MatMult(*m_basis, *m_displacements, *wksp.m_stacktmp);
  CHKERRABORT(m_comm, perr);
  wksp.scatter_stacked_to_grads();

  // build warp matrix
  std::vector<Vec*> tmps(0);
  for (auto const& vptr : wksp.m_globaltmps)
  {
    tmps.push_back(vptr.get());
  }
  Mat_unique warp = build_warp_matrix(m_comm, m_v_image_shape, image.ndim(), tmps);

  // now apply matrix to get new image data
  // first need image in natural ordering
  Vec_unique src_nat = create_unique_vec();
  perr = DMDACreateNaturalVector(*image.dmda(), src_nat.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*src_nat, "Vec_source_natural");
  perr = DMDAGlobalToNaturalBegin(*image.dmda(), *image.global_vec(), INSERT_VALUES, *src_nat);
  CHKERRABORT(m_comm, perr);
  perr = DMDAGlobalToNaturalEnd(*image.dmda(), *image.global_vec(), INSERT_VALUES, *src_nat);
  CHKERRABORT(m_comm, perr);
  // do mult
  Vec_unique tgt_nat = create_unique_vec();
  perr = DMDACreateNaturalVector(*image.dmda(), tgt_nat.get());
  CHKERRABORT(m_comm, perr);
  debug_creation(*tgt_nat, "Vec_target_natural");
  perr = MatMult(*warp, *src_nat, *tgt_nat);
  CHKERRABORT(m_comm, perr);

  // create new image and insert data in petsc ordering
  std::unique_ptr<Image> new_image = image.duplicate();
  perr =
      DMDANaturalToGlobalBegin(*image.dmda(), *tgt_nat, INSERT_VALUES, *new_image->global_vec());
  CHKERRABORT(m_comm, perr);
  perr = DMDANaturalToGlobalEnd(*image.dmda(), *tgt_nat, INSERT_VALUES, *new_image->global_vec());
  CHKERRABORT(m_comm, perr);
  return new_image;
}

std::pair<intvector, intvector> Map::get_dmda_local_extents() const
{
  initialize_dmda();

  intvector locs(3, 0);
  intvector widths(3, 0);

  PetscErrorCode perr =
      DMDAGetCorners(*map_dmda, &locs[0], &locs[1], &locs[2], &widths[0], &widths[1], &widths[2]);
  CHKERRABORT(m_comm, perr);

  return std::make_pair(locs, widths);
}

Vec_unique Map::get_dim_data_dmda_blocked(integer dim) const
{
  initialize_dmda();
  if (dim >= m_ndim)
  {
    throw std::runtime_error("Index too large for map dimensions");
  }
  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed
                    // obj.
  PetscErrorCode perr =
      DMDAGetAO(*map_dmda, &ao_petsctonat); // Destroying this would break the dmda
  CHKERRABORT(m_comm, perr);

  // Allocate temp vec
  Vec_unique tmp_data = create_unique_vec();
  DMCreateGlobalVector(*map_dmda, tmp_data.get());

  // Get extents of local data in grad array
  integer startelem, datasize;
  perr = VecGetOwnershipRange(*tmp_data, &startelem, &datasize);
  CHKERRABORT(m_comm, perr);
  datasize -= startelem;

  // Target range is then just owned part of new data vec
  IS_unique tgt_is(create_unique_is());
  perr = ISCreateStride(m_comm, datasize, startelem, 1, tgt_is.get());
  CHKERRABORT(m_comm, perr);

  // Source range is equivalent range in appropriate ordering
  IS_unique src_is(create_unique_is());
  perr = ISCreateStride(m_comm, datasize, startelem + (dim * this->size()), 1, src_is.get());
  perr = AOApplicationToPetscIS(ao_petsctonat, *src_is);
  CHKERRABORT(m_comm, perr);

  // Now create the scatter
  VecScatter_unique sct = create_unique_vecscatter();
  perr = VecScatterCreate(*m_displacements, *src_is, *tmp_data, *tgt_is, sct.get());
  CHKERRABORT(m_comm, perr);

  perr = VecScatterBegin(*sct, *m_displacements, *tmp_data, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(m_comm, perr);
  perr = VecScatterEnd(*sct, *m_displacements, *tmp_data, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(m_comm, perr);

  return tmp_data;
}

void Map::calculate_basis()
{
  // Get the full Nd basis
  floatvector scalings(m_ndim, 0.0);
  floatvector offsets(m_ndim, 0.0);

  std::transform(
      this->m_v_node_spacing.begin(), this->m_v_node_spacing.end(), scalings.begin(),
      [](floating a) -> floating { return 1 / a; });
  n_ary_transform(
      [](floating x, floating a) -> floating { return -x / a; }, offsets.begin(),
      this->m_v_offsets.begin(), this->m_v_offsets.end(), this->m_v_node_spacing.begin());
  m_basis = build_basis_matrix(
      m_comm, map_shape, m_v_image_shape, scalings, offsets, m_ndim, m_ndim + 1);

  // Now grab a 1d basis as a submatrix. Note can't do this the other way round because Petsc won't
  // allow reuse of rows/cols in MatCreateSubMatrix
  // Work out what rows, columns the rank needs to own for compatability with image and
  // displacement vectors
  /*  integer rowstart, rowsize, colstart, colsize;
    PetscErrorCode perr = VecGetOwnershipRange(*m_mask.global_vec(), &rowstart,
    &rowsize);CHKERRABORT(m_comm, perr); rowsize -= rowstart; throw std::runtime_error("need to get
    just one map length's worth of columns"); perr = MatGetOwnershipRangeColumn(*m_basis,
    &colstart, &colsize);CHKERRABORT(m_comm, perr); colsize -= colstart;

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
  PetscErrorCode perr = VecGetOwnershipRange(*m_displacements, &startrow, &endrow);
  CHKERRABORT(m_comm, perr);

  Mat_unique lapl = build_laplacian_matrix(m_comm, map_shape, startrow, endrow, m_ndim + 1);
  perr = MatTransposeMatMult(*lapl, *lapl, MAT_INITIAL_MATRIX, PETSC_DEFAULT, m_lapl.get());
  debug_creation(*m_lapl, "Mat_l_squared");
  CHKERRABORT(m_comm, perr);
}

std::pair<integer, integer> Map::get_displacement_ownershiprange() const
{
  std::pair<integer, integer> range;

  PetscErrorCode perr = VecGetOwnershipRange(*m_displacements, &range.first, &range.second);
  CHKERRABORT(m_comm, perr);

  return range;
}

const floating* Map::get_raw_data_ro() const
{
  const floating* ptr;
  PetscErrorCode perr = VecGetArrayRead(*m_displacements, &ptr);
  CHKERRABORT(m_comm, perr);
  return ptr;
}

void Map::release_raw_data_ro(const floating*& ptr) const
{
  PetscErrorCode perr = VecRestoreArrayRead(*m_displacements, &ptr);
  CHKERRABORT(m_comm, perr);
}

floatvector Map::low_corner() const
{
  floatvector corner(3, 0.);
  std::transform(m_vv_node_locs.cbegin(), m_vv_node_locs.cend(), corner.begin(),
                 [](const floatvector &v){return v.front();});
  return corner;
}
