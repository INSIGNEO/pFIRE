#include "parallelmap.hpp"
#include "serialmap.hpp"

#include <iostream>
#include <iterator>

#include "indexing.hpp"
#include "iterator_routines.hpp"

const intcoordvector MapBase::node_offset_list_2d{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {1, 1, 0}};
const intcoordvector MapBase::node_offset_list_3d{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {1, 1, 0},
                                                  {0, 0, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}};

MapBase::MapBase(const intcoord& shape, const integer& ndim, const intcoord& node_spacing, const Mask& mask,
                 const MPI_Comm& comm)
  : GridVariable(shape, ndim, MapBase::validate_map_comm(comm, mask)), _mask(mask),
    _node_spacing(node_spacing), _laplacian2(nullptr)
{
  calculate_node_locations();
}

floatcoord MapBase::coord_from_index(const intcoord& index) const
{
  return {_node_locations[0][index[0]], _node_locations[1][index[1]], _node_locations[2][index[2]]};
}

const Mat& MapBase::laplacian2() const
{
  if (!this->_laplacian2)
  {
    this->build_laplacian2();
  }

  return *this->_laplacian2;
}

void MapBase::map_print() const
{
  floating**** map_data;
  PetscErrorCode perr = DMDAVecGetArrayDOFRead(this->dmda(), this->global_vector(), &map_data);
  CHKERRXX(perr);
  for(integer idof=0; idof < this->ndof(); idof++)
  {
    std::cout << "DOF = " << idof << "\n";
    for(integer zz=0; zz < this->shape()[2]; zz++)
    {
      std::cout << "\tz = " << zz << "\n";
      for(integer yy=0; yy < this->shape()[1]; yy++)
      {
        std::cout << "\t";
        for(integer xx=0; xx < this->shape()[0]; xx++)
        {
          std::cout << map_data[zz][yy][xx][idof] << " ";
        }
        std::cout << "\n";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  perr = DMDAVecRestoreArrayDOFRead(this->dmda(), this->global_vector(), &map_data);
  CHKERRXX(perr);

}

void MapBase::build_laplacian2() const
{
  integer row_lo, row_hi;
  PetscErrorCode perr = VecGetOwnershipRange(this->displacement_vector(), &row_lo, &row_hi);
  CHKERRXX(perr);
  integer num_local_rows = row_hi - row_lo;
  integer num_global_rows = this->num_total_nodes() * this->ndof();

  auto col_offsets = calculate_von_neumann_offsets(this->shape(), this->ndof());

  intvector row_ptrs;
  intvector column_ptrs;
  floatvector data;
  integer row_ptr = 0;
  row_ptrs.push_back(row_ptr);
  for (integer curr_row = row_lo; curr_row < row_hi; curr_row++)
  {
    floating row_start = row_ptr;
    integer diag_loc = row_ptr;
    for (const auto& offset : col_offsets)
    {
      integer column_ptr = curr_row + offset;
      if (column_ptr < 0 || column_ptr >= num_global_rows)
      {
        continue;
      }
      column_ptrs.push_back(column_ptr);
      data.push_back(1.);
      if (offset == 0)
      {
        diag_loc = row_ptr;
      }
      row_ptr++;
    }
    data[diag_loc] = 1 - (row_ptr - row_start);
    row_ptrs.push_back(row_ptr);
  }

  MPI_Comm vec_comm;
  PetscObjectGetComm(reinterpret_cast<PetscObject>(this->displacement_vector()), &vec_comm);
  Mat_unique lapl = create_unique_mat();
  perr = MatCreateMPIAIJWithArrays(vec_comm, num_local_rows, num_local_rows, // Square subsection
                                   num_global_rows, num_global_rows,             // cross-check
                                   row_ptrs.data(), column_ptrs.data(),          // csr data
                                   data.data(), lapl.get());
  CHKERRXX(perr);

  _laplacian2 = create_unique_mat();
  perr = MatTransposeMatMult(*lapl, *lapl, MAT_INITIAL_MATRIX, PETSC_DEFAULT, _laplacian2.get());
  CHKERRXX(perr);
}

intcoord MapBase::calculate_map_shape(const intcoord& mask_shape, const intcoord& spacing)
{
  intcoord map_shape = {1, 1, 1};

  for (size_t idx = 0; idx < map_shape.size(); idx++)
  {
    if (mask_shape[idx] == 1)
    {
      continue;
    }
    auto n_nodes = static_cast<integer>(std::ceil((mask_shape[idx] + 1.) / spacing[idx]));
    map_shape[idx] = (n_nodes > 3) ? n_nodes : 3;
  }
  return map_shape;
}

std::pair<intcoord, intcoord> MapBase::get_pixel_neighbourhood(const intcoord& map_node) const
{
  intcoord cimlo = {0, 0, 0};
  intcoord cimhi = {0, 0, 1};
  for (integer idim = 0; idim < this->ndim(); idim++)
  {
    floating ctr = this->node_locs()[idim][map_node[idim]];
    cimlo[idim] = std::ceil(ctr - this->spacing()[idim]);
    cimhi[idim] = std::floor(ctr + this->spacing()[idim]);
  }
  clamp_location_lo(cimlo, this->mask().index_min());
  clamp_location_hi(cimhi, this->mask().shape());

  return std::make_pair(cimlo, cimhi);
}

std::unique_ptr<MapBase> MapBase::make_map_for_mask(const Mask& mask, const intcoord& node_spacing)
{
  // Determine shape from mask and nodespacing
  intcoord map_shape = MapBase::calculate_map_shape(mask.shape(), node_spacing);

  integer ndof = mask.ndim() + 1;

  // Calculate what parallel shape *would* be
  intcoord rank_distribution = GridVariable::find_rank_distribution(mask.shape(), mask.commsize());

  std::unique_ptr<MapBase> map;
  if (all_true(map_shape.cbegin(), map_shape.cbegin() + mask.ndim(), rank_distribution.cbegin(),
               rank_distribution.cend() + mask.ndim(),
               [](integer m, integer n) -> bool { return m > (2 * n + 1); }))
  {
    map = std::make_unique<ParallelMap>(map_shape, ndof, node_spacing, mask);
  }
  else
  {
    map = std::make_unique<SerialMap>(map_shape, ndof, node_spacing, mask);
  }

  return map;
}

MPI_Comm MapBase::validate_map_comm(const MPI_Comm& comm, const Mask& mask)
{
  if (comm == MPI_COMM_SELF || comm == mask.comm())
  {
    return comm;
  }

  throw InternalError("Attempted to create map with invalid MPI_Comm");
}

void MapBase::calculate_node_locations()
{
  floatcoord map_width = {0., 0., 0.};
  n_ary_transform([](floating n_nodes, floating spacing) -> floating { return spacing * (n_nodes - 1); },
                  map_width.begin(), this->shape().begin(), this->shape().end(), this->spacing().begin());

  // if (any_true(map_width.cbegin(), map_width.cend(), this->mask().shape().cbegin(),
  // std::less<>()))
  //{
  // throw InternalError("Map too small for associated mask!");
  //}

  const floating half = 0.5;
  floatcoord margin = half * (map_width - this->mask().shape());
  floatcoord lo_corner = (this->mask().index_min() - margin) - half; // image is cell centered but the
                                                                        // map is node centered

  _node_locations.reserve(this->shape().size());
  for (size_t idim = 0; idim < this->shape().size(); idim++)
  {
    floatvector nodes;
    if (this->shape()[idim] == 1)
    {
      nodes.push_back(this->mask().index_min()[idim]);
    }
    else
    {
      floating curr_loc = lo_corner[idim];
      for (integer iidx = 0; iidx < this->shape()[idim]; iidx++)
      {
        nodes.push_back(curr_loc);
        curr_loc += this->spacing()[idim];
      }
    }
    _node_locations.push_back(nodes);
  }
}

floatcoordvector MapBase::map_local_coordinates(const intcoordvector& loclist) const
{
  const intcoordvector* offsets;
  if (this->ndim() == 2)
  {
    offsets = &MapBase::node_offset_list_2d;
  }
  else
  {
    offsets = &MapBase::node_offset_list_3d;
  }

  // Acquire the relevant vectors...
  floating**** map_data;
  PetscErrorCode perr = DMDAVecGetArrayDOFRead(this->dmda(), this->local_vector(), &map_data);
  CHKERRXX(perr);

  // For each node in list, locate neighbours in previous and add interpolation
  floatcoordvector mapped_locs;
  mapped_locs.reserve(loclist.size());
  for (const auto& curr_coord : loclist)
  {
    floatcoord mapped_loc = curr_coord;
    intcoord corner_node = this->get_cell_corner(curr_coord);
    for (const auto& offset : *offsets)
    {
      intcoord map_node = corner_node + offset;
      /*if (map_node >= previous_map.shape())
      {
        continue;
      }*/
      floatcoord map_coord = this->coord_from_index(map_node);

      // Interpolation coefficients are given by
      // \prod_{i=x,y,z} 1 - | a_i - b_i|/s_i where a_i are the previous node coordinates
      // b_i are the new node coordinates and s_i is the node spacing of the *previous* map
      floating interp_coeff = calculate_scaled_basis_coefficient(
          curr_coord.cbegin(), curr_coord.cend(), map_coord.cbegin(), this->spacing().cbegin());
      // recall idof[0] is intensity so to update dimensional coordinates start from 1
      for (integer idof = 0; idof < this->ndof(); idof++)
      {
        mapped_loc[idof] -= interp_coeff * map_data[map_node[2]][map_node[1]][map_node[0]][idof+1];
      }
    }
    clamp_location(mapped_loc, floatcoord({0,0,0}), this->mask().shape() - floatcoord({1,1,1}));
    mapped_locs.push_back(mapped_loc);
  }

  perr = DMDAVecRestoreArrayDOFRead(this->dmda(), this->local_vector(), &map_data);
  CHKERRXX(perr);

  return mapped_locs;
}
