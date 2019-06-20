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

#ifndef MAPMANAGER_HPP
#define MAPMANAGER_HPP

#include "types.hpp"
#include "mapbase.hpp"

class MapManager
{
public:
  MapManager(const intcoord& target_spacing, const Mask& mask);

  integer step_number() const {return _step_number;}

  //! Given a target spacing and mask, calculate all spacings needed in reverse order
  static intcoordvector calculate_spacing_list(const intcoord& mask_shape, const intcoord& target_spacing);

  //! Get the current map
  std::shared_ptr<MapBase> get_current_map() const { return _current_map;};

  //! Get the current map
  const intcoord& get_current_spacing() const { return _spacing_list.back();};

  //! Get the next most refined map by interpolating the current map
  std::shared_ptr<MapBase> get_next_map();

private:


  const Mask& _mask;
  intcoordvector _spacing_list;  // _spacing_list.back() should at all times refer to the *current* spacing
  std::shared_ptr<MapBase> _current_map;
  integer _step_number;

};

#endif //MAPMANAGER_HPP
