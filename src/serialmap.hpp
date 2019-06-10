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

#ifndef SERIALMAP_HPP
#define SERIALMAP_HPP

#include "mapbase.hpp"

class SerialMap: public MapBase
{
public:
  SerialMap(const intcoord &shape, const integer& ndof, const intcoord& node_spacing, const Mask& mask)
    : MapBase(shape, ndof, node_spacing, mask, MPI_COMM_SELF), _map_shared_vec(create_unique_vec())
  {
    _shared_comm = mask.comm();
    PetscErrorCode perr = VecCreateMPI(_shared_comm, PETSC_DECIDE, this->num_total_nodes() * this->ndof(),
                                       _map_shared_vec.get());
    CHKERRXX(perr);
    
    perr = VecSet(*_map_shared_vec, 0.);
    CHKERRXX(perr);
  }

  //! Copy vector data in from compatible vector
  void copy_data_from(const SerialMap& src);

  bool parallel_layout() const override {return false;}

  const Vec& displacement_vector() const override {return *_map_shared_vec;}

  void interpolate_from(const MapBase& previous_map) override; 

  void update(const Vec& delta) override;

private:

  void interpolate_serial_to_serial(const SerialMap& previous_map);
  void interpolate_parallel_to_serial(const ParallelMap& previous_map);

  Vec_unique _map_shared_vec;

  MPI_Comm _shared_comm;

  friend class ParallelMap;
};

#endif //SERIALMAP_HPP
