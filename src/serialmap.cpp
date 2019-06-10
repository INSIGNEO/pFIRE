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

#include "serialmap.hpp"

#include "exceptions.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "parallelmap.hpp"

void SerialMap::interpolate_from(const MapBase& previous_map)
{
  if (previous_map.parallel_layout())
  {
    this->interpolate_parallel_to_serial(dynamic_cast<const ParallelMap&>(previous_map));
  }
  else
  {
    this->interpolate_serial_to_serial(dynamic_cast<const SerialMap&>(previous_map));
  }
}

void SerialMap::interpolate_serial_to_serial(const SerialMap& previous_map)
{
  PetscPrintf(this->comm(), "Interpolating map, serial->serial\n");
  // All data is local, so no MPI complications here.

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

  for (zz = 0; zz < this->shape()[2]; zz++)
  {
    for (yy = 0; yy < this->shape()[1]; yy++)
    {
      for (xx = 0; xx < this->shape()[0]; xx++)
      {
        floatcoord curr_coord = this->coord_from_index(curr_node);
        intcoord corner_node = previous_map.get_cell_corner(curr_coord);
        for (const auto& offset : *offsets)
        {
          intcoord prev_node = corner_node + offset;
          /*if (prev_node >= previous_map.shape())
          {
            continue;
          }*/
          floatcoord previous_coord = previous_map.coord_from_index(prev_node);

          // Interpolation coefficients are given by
          // \prod_{i=x,y,z} 1 - | a_i - b_i|/s_i where a_i are the previous node coordinates
          // b_i are the new node coordinates and s_i is the node spacing of the *previous* map
          floating interp_coeff = calculate_scaled_basis_coefficient(curr_coord.cbegin(), curr_coord.cend(),
                                                                     previous_coord.cbegin(),
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

void SerialMap::interpolate_parallel_to_serial(const ParallelMap& previous_map __attribute((unused)))
{
  // If we are trying to go from parallel to serial then we are trying to coarsen
  // which is not needed in the general flow of pFIRE
  throw InternalError("Coarsening of maps is not implemented");
}

void SerialMap::copy_data_from(const SerialMap& src)
{
  if (this->dmda() != src.dmda())
  {
    throw InternalError("GridVariables must share a DMDA", __FILE__, __LINE__);
  }

  // Do the vector copy
  PetscErrorCode perr = VecCopy(src.global_vector(), this->global_vector());
  CHKERRXX(perr);
  perr = VecCopy(src.displacement_vector(), this->displacement_vector());
  CHKERRXX(perr);

  this->update_local_vector();
}

void SerialMap::update(const Vec& delta)
{
  PetscErrorCode perr = VecAXPY(this->displacement_vector(), 1., delta);
  CHKERRXX(perr);

  integer owned_lo, owned_hi;
  perr = VecGetOwnershipRange(this->displacement_vector(), &owned_lo, &owned_hi);
  CHKERRXX(perr);

  perr = VecSet(this->global_vector(), 0.);

  floating**** map_data;
  perr = DMDAVecGetArrayDOF(this->dmda(), this->global_vector(), &map_data);
  CHKERRXX(perr);
  const floating* vec_data;
  perr = VecGetArrayRead(this->displacement_vector(), &vec_data);
  CHKERRXX(perr);

  for (integer iidx = owned_lo; iidx < owned_hi; iidx++)
  {
    integer idof = iidx % this->ndof();
    intcoord iloc = unravel(iidx / this->ndof(), this->shape());
    map_data[iloc[2]][iloc[1]][iloc[0]][idof] += vec_data[iidx-owned_lo];
  }

  perr = VecRestoreArrayRead(this->displacement_vector(), &vec_data);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArrayDOF(this->dmda(), this->global_vector(), &map_data);
  CHKERRXX(perr);

  floating* vec_data_rw;
  perr = VecGetArray(this->global_vector(), &vec_data_rw);
  CHKERRXX(perr);
  MPI_Allreduce(MPI_IN_PLACE, vec_data_rw, this->ndof() * this->num_total_nodes(), MPIU_SCALAR, MPI_SUM,
                this->_shared_comm);
  perr = VecRestoreArray(this->global_vector(), &vec_data_rw);
  CHKERRXX(perr);


  this->update_local_vector();
}
