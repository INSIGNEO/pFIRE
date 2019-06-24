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

#ifndef PARALLELMAP_HPP
#define PARALLELMAP_HPP

#include "mapbase.hpp"

class ParallelMap: public MapBase {
public:
  ParallelMap(const intcoord& shape, const integer& ndof, const intcoord& node_spacing, const Mask& mask)
    : MapBase(shape, ndof, node_spacing, mask, mask.comm())
  {
  }

  bool parallel_layout() const override { return true; }

  const Vec& displacement_vector() const { return this->global_vector(); }

  void interpolate_from(const MapBase& previous_map) override;

  void update(const Vec& delta) override;

private:
  void interpolate_serial_to_parallel(const SerialMap& previous_map);
  void interpolate_parallel_to_parallel(const ParallelMap& previous_map);

  friend class SerialMap;
};

#endif // PARALLELMAP_HPP
