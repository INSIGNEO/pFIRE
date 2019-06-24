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

#include "parallelmap.hpp"

#include "serialmap.hpp"
#include "exceptions.hpp"
#include "indexing.hpp"

void ParallelMap::interpolate_from(const MapBase& previous_map)
{
  if(previous_map.parallel_layout())
  {
    this->interpolate_parallel_to_parallel(dynamic_cast<const ParallelMap&>(previous_map));
  }
  else
  {
    this->interpolate_serial_to_parallel(dynamic_cast<const SerialMap&>(previous_map));
  }
}

void ParallelMap::interpolate_serial_to_parallel(const SerialMap& previous_map)
{
  // All source data is local, so no MPI comm needed, just act on locally owned target data
  auto local_range = this->owned_range();
  intcoord& local_lo = local_range.first;
  intcoord& local_hi = local_range.second;

  // Zero vector data first...
  PetscErrorCode perr = VecSet(this->global_vector(), 0.);

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
  floating**** current_data;
  perr = DMDAVecGetArrayDOF(this->dmda(), this->global_vector(), &current_data);
  CHKERRXX(perr);
  floating**** previous_data;
  perr = DMDAVecGetArrayDOF(previous_map.dmda(), previous_map.global_vector(), &previous_data);
  CHKERRXX(perr);

  // For each node in current, locate neighbours in previous and add interpolation
  intcoord curr_node;
  integer& xx = curr_node[0];
  integer& yy = curr_node[1];
  integer& zz = curr_node[2];

  for (zz = local_lo[2]; zz < local_hi[2]; zz++)
  {
    for (yy = local_lo[1]; yy < local_hi[1]; yy++)
    {
      for (xx = local_lo[0]; xx < local_hi[0]; xx++)
      {
        floatcoord curr_coord = this->coord_from_index(curr_node);
        intcoord corner_node = previous_map.get_cell_corner(curr_coord);
        for (const auto& offset : *offsets)
        {
          intcoord prev_node = corner_node + offset;
          floatcoord previous_coord = previous_map.coord_from_index(prev_node);

          // Interpolation coefficients are given by
          // \prod_{i=x,y,z} 1 - | a_i - b_i|/s_i where a_i are the previous node coordinates
          // b_i are the new node coordinates and s_i is the node spacing of the *previous* map
          floating interp_coeff = calculate_scaled_basis_coefficient(
              curr_coord.cbegin(), curr_coord.cend(), previous_coord.cbegin(),
              previous_map.spacing().cbegin());
          for (integer idof = 0; idof < this->ndof(); idof++)
          {
            current_data[zz][yy][xx][idof] +=
                interp_coeff * previous_data[prev_node[2]][prev_node[1]][prev_node[0]][idof];
          }
        }
      }
    }
  }
  perr = DMDAVecRestoreArrayDOF(this->dmda(), this->global_vector(), &current_data);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArrayDOF(previous_map.dmda(), previous_map.global_vector(), &previous_data);
  CHKERRXX(perr);
  this->update_local_vector();
}


void ParallelMap::interpolate_parallel_to_parallel(const ParallelMap& previous_map)
{
  // In general this needs a horrible all-to-all to guarantee we get all the data, but in practice 
  // provided we are refining resolution, then we will at most need the first neighbour from other ranks
  // which is taken care of for us by the DMDA, so we just make sure the local data is up to date on the
  // previous map and use that to update the global data on the new map
  
  // Guard against attempting to coarsen
  if(this->num_total_nodes() < previous_map.num_total_nodes())
  {
    throw InternalError("Coarsening of maps is not implemented");
  }

  // Ensure local data is up to date
  previous_map.update_local_vector();

  auto local_range = this->owned_range();
  intcoord& local_lo = local_range.first;
  intcoord& local_hi = local_range.second;

  // Zero vector data first...
  PetscErrorCode perr = VecSet(this->global_vector(), 0.);

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
  floating**** current_data;
  perr = DMDAVecGetArrayDOF(this->dmda(), this->global_vector(), &current_data);
  CHKERRXX(perr);
  floating**** previous_data;
  perr = DMDAVecGetArrayDOF(previous_map.dmda(), previous_map.local_vector(), &previous_data);
  CHKERRXX(perr);

  // For each node in current, locate neighbours in previous and add interpolation
  intcoord curr_node;
  integer& xx = curr_node[0];
  integer& yy = curr_node[1];
  integer& zz = curr_node[2];

  for (zz = local_lo[2]; zz < local_hi[2]; zz++)
  {
    for (yy = local_lo[1]; yy < local_hi[1]; yy++)
    {
      for (xx = local_lo[0]; xx < local_hi[0]; xx++)
      {
        floatcoord curr_coord = this->coord_from_index(curr_node);
        intcoord corner_node = previous_map.get_cell_corner(curr_coord);
        for (const auto& offset : *offsets)
        {
          intcoord prev_node = corner_node + offset;
          floatcoord previous_coord = previous_map.coord_from_index(prev_node);

          // Interpolation coefficients are given by
          // \prod_{i=x,y,z} 1 - | a_i - b_i|/s_i where a_i are the previous node coordinates
          // b_i are the new node coordinates and s_i is the node spacing of the *previous* map
          floating interp_coeff = calculate_scaled_basis_coefficient(
              curr_coord.cbegin(), curr_coord.cend(), previous_coord.cbegin(),
              previous_map.spacing().cbegin());
          for (integer idof = 0; idof < this->ndof(); idof++)
          {
            current_data[zz][yy][xx][idof] +=
                interp_coeff * previous_data[prev_node[2]][prev_node[1]][prev_node[0]][idof];
          }
        }
      }
    }
  }
  perr = DMDAVecRestoreArrayDOF(this->dmda(), this->global_vector(), &current_data);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArrayDOF(previous_map.dmda(), previous_map.local_vector(), &previous_data);
  CHKERRXX(perr);

  this->update_local_vector();
}

void ParallelMap::update(const Vec& delta)
{
  PetscErrorCode perr = VecAXPY(this->displacement_vector(), 1., delta);
  CHKERRXX(perr);

  this->update_local_vector();
}
