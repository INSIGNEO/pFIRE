//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include "map.hpp"

#include <petscis.h>
#include <petscvec.h>

#include "basis.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "laplacian.hpp"
#include "math_utils.hpp"
#include "petsc_helpers.hpp"
#include "workspace.hpp"
#include "math_utils.hpp"

#include "iterator_routines.hpp"

Map::Map(floatvector node_spacing, const Mask& mask)
  : m_comm(mask.comm()), m_mask(mask), m_ndim(mask.ndim()),
    m_v_node_spacing(std::move(node_spacing)), m_v_image_shape(mask.shape()),
    m_basis(create_unique_mat()), m_lapl(create_unique_mat()),
    m_displacements(create_unique_vec()), map_dmda(create_unique_dm())
{
  calculate_node_info();
  calculate_basis();
  initialize_dmda();

  apply_mask_to_basis();

  // initialize displacement storage
  alloc_displacements();
  calculate_laplacian();
}


void Map::apply_mask_to_basis()
{
  // Need to apply mask to both image and map "side" of basis
  // First stack basis to depth of map
  Vec_unique stacked_basis = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(*m_basis, nullptr, stacked_basis.get());
  CHKERRXX(perr);

  repeat_stack(*m_mask.global_vec(), *stacked_basis);

  // For map side determine from mask what nodes are relevant, mask out others
  Vec_unique map_side = calculate_map_mask(*stacked_basis);

  // Use matdiagonalscale and apply both at once
  // Mask applies directly to image side of basis

  perr = MatDiagonalScale(*m_basis, *stacked_basis, *map_side);
  CHKERRXX(perr);
}

Vec_unique Map::calculate_map_mask(Vec& stacked_mask )
{
  // Create map mask from image mask and basis
  Vec_unique map_mask = create_unique_vec();
  PetscErrorCode perr = MatCreateVecs(*m_basis, map_mask.get(), nullptr);
  CHKERRXX(perr);
  perr = MatMultTranspose(*m_basis, stacked_mask, *map_mask);
  CHKERRXX(perr);

  
  // Dilate to include surrounding nodes
  // First need 1 copy in global layout
  Vec_unique map_mask_global = create_unique_vec();
  perr = DMCreateGlobalVector(*map_dmda, map_mask_global.get());
  CHKERRXX(perr);
  copy_nth_from_stack_nat_to_petsc(*map_mask_global, *map_mask, *map_dmda, 0);

  // Now dilate
  dilate_dmda(*map_dmda, *map_mask_global);
  binarize_vector(*map_mask_global, 0.01);

  // Now map back:
  repeat_stack_petsc_to_nat(*map_mask_global, *map_mask, *map_dmda);

  return map_mask;
}

void Map::update(const Vec& delta_vec)
{
  PetscErrorCode perr = VecAXPY(*m_displacements, 1, delta_vec);
  CHKERRXX(perr);
}

std::unique_ptr<Map> Map::interpolate(const floatvector& new_spacing)
{
  std::unique_ptr<Map> new_map(new Map(new_spacing, this->m_mask));

  floatvector scalings(m_ndim, 0.0);
  floatvector offsets(m_ndim, 0.0);

  n_ary_transform(std::divides<>(), scalings.begin(), new_map->m_v_node_spacing.begin(),
      new_map->m_v_node_spacing.end(), this->m_v_node_spacing.begin());
  n_ary_transform([](floating x, floating y, floating a) -> floating { return (x - y) / a; },
      offsets.begin(), new_map->m_v_offsets.begin(), new_map->m_v_offsets.end(),
      this->m_v_offsets.begin(), this->m_v_node_spacing.begin());

  Mat_unique interp = build_basis_matrix(
      m_comm, map_shape, new_map->map_shape, scalings, offsets, m_ndim, m_ndim + 1);

  PetscErrorCode perr = MatMult(*interp, *m_displacements, *new_map->m_displacements);
  CHKERRXX(perr);

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
  perr = DMDACreate3d(m_comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, // BCs
      DMDA_STENCIL_STAR,                        // stencil shape
      map_shape[0], map_shape[1], map_shape[2], // global mesh shape
      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, // ranks per dim
      dof_per_node, stencil_width,              // dof per node, stencil size
      nullptr, nullptr, nullptr,                // partition sizes nullptr -> petsc chooses
      map_dmda.get());
  CHKERRXX(perr);

  perr = DMSetUp(*(map_dmda));
  CHKERRXX(perr);
}

void Map::alloc_displacements() { MatCreateVecs(*m_basis, m_displacements.get(), nullptr); }

intvector Map::calculate_map_shape(intvector const& image_shape, floatvector const& nodespacing)
{
  if (image_shape.size() != nodespacing.size())
  {
    throw InternalError("Image and nodespacing dimensions must match", __FILE__, __LINE__);
  }

  // want always to have odd number of nodes so find num nodes for each half,
  // multiply by two and subtract one to get total nodes
  // N.B number of nodes = 1 + number of spaces, so:
  intvector nnod(image_shape.size());
  std::transform(image_shape.cbegin(), image_shape.cend(), nodespacing.cbegin(), nnod.begin(),
      [](integer is, floating ns) -> integer {
        return static_cast<integer>(1 + 2 * std::ceil((is / 2. - 0.5) / ns));
      });
  return nnod;
}

void Map::calculate_node_info()
{
  map_shape = calculate_map_shape(m_v_image_shape, m_v_node_spacing);
  // calculate self size and offset
  for (uinteger idim = 0; idim < m_ndim; idim++)
  {
    floating lo = 0.5 * (m_v_image_shape[idim] - (map_shape[idim] - 1) * m_v_node_spacing[idim]);
    floatvector nodes(map_shape[idim], 0.0);
    std::generate(nodes.begin(), nodes.end(), [x = lo, y = m_v_node_spacing[idim]]() mutable {
      x += y;
      return x - y;
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
  CHKERRXX(perr);
  wksp.scatter_stacked_to_grads_noreorder();

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
  CHKERRXX(perr);
  debug_creation(*src_nat, "Vec_source_natural");
  perr = DMDAGlobalToNaturalBegin(*image.dmda(), *image.global_vec(), INSERT_VALUES, *src_nat);
  CHKERRXX(perr);
  perr = DMDAGlobalToNaturalEnd(*image.dmda(), *image.global_vec(), INSERT_VALUES, *src_nat);
  CHKERRXX(perr);
  // do mult
  Vec_unique tgt_nat = create_unique_vec();
  perr = DMDACreateNaturalVector(*image.dmda(), tgt_nat.get());
  CHKERRXX(perr);
  debug_creation(*tgt_nat, "Vec_target_natural");
  perr = MatMult(*warp, *src_nat, *tgt_nat);
  CHKERRXX(perr);

  // create new image and insert data in petsc ordering
  std::unique_ptr<Image> new_image = Image::duplicate(image);
  perr =
      DMDANaturalToGlobalBegin(*image.dmda(), *tgt_nat, INSERT_VALUES, *new_image->global_vec());
  CHKERRXX(perr);
  perr = DMDANaturalToGlobalEnd(*image.dmda(), *tgt_nat, INSERT_VALUES, *new_image->global_vec());
  CHKERRXX(perr);
  return new_image;
}

std::pair<intvector, intvector> Map::get_dmda_local_extents() const
{
  initialize_dmda();

  intvector locs(3, 0);
  intvector widths(3, 0);
  PetscErrorCode perr =
      DMDAGetCorners(*map_dmda, &locs[0], &locs[1], &locs[2], &widths[0], &widths[1], &widths[2]);
  CHKERRXX(perr);

  return std::make_pair(locs, widths);
}

Vec_unique Map::get_raw_data_row_major(uinteger dim) const
{
  // scatter relevant dim to dmda vec, reordering natural to global
  Vec_unique tmp_vec(get_single_dim_petsc(dim));

  // get DMDA local sizes
  intvector locs(3, 0);
  intvector widths(3, 0);
  PetscErrorCode perr =
      DMDAGetCorners(*map_dmda, &locs[0], &locs[1], &locs[2], &widths[0], &widths[1], &widths[2]);
  CHKERRXX(perr);
  integer startelem, localsize;
  perr = VecGetOwnershipRange(*tmp_vec, &startelem, &localsize);
  localsize -= startelem;

  // source is simple stride
  IS_unique src_is = create_unique_is();
  perr = ISCreateStride(m_comm, localsize, startelem, 1, src_is.get());
  CHKERRXX(perr);

  // do rowmaj -> colmaj scatter as in Image
  intvector cmidxn(localsize);
  std::iota(cmidxn.begin(), cmidxn.end(), 0);
  std::transform(
      cmidxn.begin(), cmidxn.end(), cmidxn.begin(), [widths, startelem](integer idx) -> integer {
        return idx_cmaj_to_rmaj(idx, widths) + startelem;
      });
  IS_unique tgt_is = create_unique_is();
  perr = ISCreateGeneral(m_comm, localsize, cmidxn.data(), PETSC_USE_POINTER, tgt_is.get());
  CHKERRXX(perr);

  // now do scatter
  Vec_unique cmaj_vec(create_unique_vec());
  perr = VecDuplicate(*tmp_vec, cmaj_vec.get());
  CHKERRXX(perr);
  VecScatter_unique sct(create_unique_vecscatter());
  perr = VecScatterCreate(*tmp_vec, *src_is, *cmaj_vec, *tgt_is, sct.get());
  CHKERRXX(perr);
  perr = VecScatterBegin(*sct, *tmp_vec, *cmaj_vec, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);
  perr = VecScatterEnd(*sct, *tmp_vec, *cmaj_vec, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);

  return cmaj_vec;
}

/*
Vec_unique Map::get_single_dim_petsc(uinteger dim) const
{
  initialize_dmda();
  if (dim >= m_ndim)
  {
    throw InternalError("Index too large for map dimensions", __FILE__, __LINE__);
  }
  // Allocate temp vec
  Vec_unique tmp_data = create_unique_vec();
  PetscErrorCode perr = DMCreateGlobalVector(*map_dmda, tmp_data.get());
  CHKERRXX(perr);

  AO ao_petsctonat; // N.B this is not going to be a leak, we are just borrowing a Petsc managed
                    // obj.
  perr = DMDAGetAO(*map_dmda, &ao_petsctonat); // Destroying this would break the dm
  CHKERRXX(perr);

  // Target range is owned part of new vec
  integer startelem, localsize;
  perr = VecGetOwnershipRange(*tmp_data, &startelem, &localsize);
  localsize -= startelem;
  IS_unique tgt_is(create_unique_is());
  perr = ISCreateStride(m_comm, localsize, startelem, 1, tgt_is.get());
  CHKERRXX(perr);
  perr = AOPetscToApplicationIS(ao_petsctonat, *tgt_is);
  CHKERRXX(perr);

<<<<<<< HEAD
  // Source range is equivalent range offset to dimension
=======
  // Source range is owned part of vec offset by dimension length
  integer offset = dim * this->size();
>>>>>>> develop
  IS_unique src_is(create_unique_is());
  perr = ISCreateStride(m_comm, localsize, startelem + offset, 1, src_is.get());
  CHKERRXX(perr);

  VecScatter_unique sct = create_unique_vecscatter();
  perr = VecScatterCreate(*m_displacements, *src_is, *tmp_data, *tgt_is, sct.get());
  CHKERRXX(perr);

  perr = VecScatterBegin(*sct, *m_displacements, *tmp_data, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);
  perr = VecScatterEnd(*sct, *m_displacements, *tmp_data, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);

  return tmp_data;
}*/

Vec_unique Map::get_single_dim_petsc(uinteger dim) const
{
  initialize_dmda();
  if (dim >= m_ndim)
  {
    throw InternalError("Index too large for map dimensions", __FILE__, __LINE__);
  }
  // Get dim data
  Vec_unique nat_data = get_single_dim_natural(dim);
  // Allocate temp vec
  Vec_unique tmp_data(create_unique_vec());
  PetscErrorCode perr = DMCreateGlobalVector(*map_dmda, tmp_data.get());
  CHKERRXX(perr);

  DMDANaturalToGlobalBegin(*map_dmda, *nat_data, INSERT_VALUES, *tmp_data);
  DMDANaturalToGlobalEnd(*map_dmda, *nat_data, INSERT_VALUES, *tmp_data);

  return tmp_data;
}

Vec_unique Map::get_single_dim_natural(uinteger dim) const
{
  initialize_dmda();
  if (dim >= m_ndim)
  {
    throw InternalError("Index too large for map dimensions", __FILE__, __LINE__);
  }
  // Allocate temp vec
  Vec_unique tmp_data = create_unique_vec();
  PetscErrorCode perr = DMDACreateNaturalVector(*map_dmda, tmp_data.get());
  CHKERRXX(perr);

  // Target range is owned part of new vec
  integer startelem, localsize;
  perr = VecGetOwnershipRange(*tmp_data, &startelem, &localsize);
  localsize -= startelem;
  IS_unique tgt_is(create_unique_is());
  perr = ISCreateStride(m_comm, localsize, startelem, 1, tgt_is.get());
  CHKERRXX(perr);

  // Source range is owned part of vec offset by dimension length
  integer offset = dim * this->size();
  IS_unique src_is(create_unique_is());
  perr = ISCreateStride(m_comm, localsize, startelem + offset, 1, src_is.get());

  VecScatter_unique sct = create_unique_vecscatter();
  perr = VecScatterCreate(*m_displacements, *src_is, *tmp_data, *tgt_is, sct.get());
  CHKERRXX(perr);

  perr = VecScatterBegin(*sct, *m_displacements, *tmp_data, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);
  perr = VecScatterEnd(*sct, *m_displacements, *tmp_data, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRXX(perr);

  return tmp_data;
}

void Map::calculate_basis()
{
  // Get the full Nd basis
  floatvector scalings(this->m_v_node_spacing.size(), 0.0);
  floatvector offsets(this->m_v_node_spacing.size(), 0.0);

  std::transform(this->m_v_node_spacing.begin(), this->m_v_node_spacing.end(), scalings.begin(),
      [](floating a) -> floating { return 1 / a; });
  n_ary_transform([](floating x, floating a) -> floating { return -x / a; }, offsets.begin(),
      this->m_v_offsets.begin(), this->m_v_offsets.end(), this->m_v_node_spacing.begin());
  m_basis = build_basis_matrix(
      m_comm, map_shape, m_v_image_shape, scalings, offsets, m_ndim, m_ndim + 1);
}

void Map::calculate_laplacian()
{
  integer startrow, endrow;
  PetscErrorCode perr = VecGetOwnershipRange(*m_displacements, &startrow, &endrow);
  CHKERRXX(perr);

  Mat_unique lapl = build_laplacian_matrix(m_comm, map_shape, startrow, endrow, m_ndim + 1);
  perr = MatTransposeMatMult(*lapl, *lapl, MAT_INITIAL_MATRIX, PETSC_DEFAULT, m_lapl.get());
  debug_creation(*m_lapl, "Mat_l_squared");
  CHKERRXX(perr);
}

std::pair<integer, integer> Map::get_displacement_ownershiprange() const
{
  std::pair<integer, integer> range;

  PetscErrorCode perr = VecGetOwnershipRange(*m_displacements, &range.first, &range.second);
  CHKERRXX(perr);

  return range;
}

const floating* Map::get_raw_data_ro() const
{
  const floating* ptr;
  PetscErrorCode perr = VecGetArrayRead(*m_displacements, &ptr);
  CHKERRXX(perr);
  return ptr;
}

void Map::release_raw_data_ro(const floating*& ptr) const
{
  PetscErrorCode perr = VecRestoreArrayRead(*m_displacements, &ptr);
  CHKERRXX(perr);
}

floatvector Map::low_corner() const
{
  floatvector corner(3, 0.);
  std::transform(m_vv_node_locs.cbegin(), m_vv_node_locs.cend(), corner.begin(),
      [](const floatvector& v) { return v.front(); });
  return corner;
}

void Map::print_dm_info()
{
  int rank;
  MPI_Comm_rank(m_comm, &rank);
  auto corners = this->get_dmda_local_extents();
  PetscSynchronizedPrintf(m_comm, "Rank %i: corner: %i %i %i: extent: %i %i %i\n", rank,
      corners.first[0], corners.first[1], corners.first[2], corners.second[0], corners.second[1],
      corners.second[2]);
  PetscSynchronizedFlush(m_comm, PETSC_STDOUT);
}
