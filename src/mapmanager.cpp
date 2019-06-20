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

#include "mapmanager.hpp"

#include "iterator_routines.hpp"

MapManager::MapManager(const intcoord& target_spacing, const Mask& mask)
  : _mask(mask), _spacing_list(calculate_spacing_list(mask.shape(), target_spacing)),
    _current_map(MapBase::make_map_for_mask(_mask, _spacing_list.back())),
    _step_number(0)
{
  PetscPrintf(mask.comm(), "Node spacing sequence: \n");
  for (auto& spacing: reverse_iterator_adapter(_spacing_list))
  {
    PetscPrintf(mask.comm(), "\t%d, %d, %d,\n", spacing[0], spacing[1], spacing[2]);
  }
}

intcoordvector MapManager::calculate_spacing_list(const intcoord& mask_shape, const intcoord& target_spacing)
{
  intcoord spacing = target_spacing;
  intcoordvector spacing_list;
  spacing_list.push_back(spacing);
  while(any_true(spacing.cbegin(), spacing.cend(), mask_shape.cbegin(),
                 [](const integer& spc, const integer& width) -> bool {return width > 2*spc + 1;}))
  {
    for(auto& space: spacing)
    {
      space *= 2;
    }
    spacing_list.push_back(spacing);
  }
  return spacing_list;
}


//! Note that spacing_list.back() should at all times refer to the current spacing,
//  so pop after creating new map - notionally more safe if exception thrown while creating new map
//  and allows us to return nullptr when incrementing off the end while get_current_map() still returns the
//  final valid map
std::shared_ptr<MapBase> MapManager::get_next_map(){
  if(_spacing_list.crbegin()+1 == _spacing_list.crend())
  {
    return nullptr;
  }

  // create the new map and interpolate data from previous
  std::shared_ptr<MapBase> new_map = MapBase::make_map_for_mask(_mask, _spacing_list.crbegin()[1]);
  new_map->interpolate_from(*_current_map);

  _step_number++;
  _current_map = new_map;
  _spacing_list.pop_back();
  return new_map;
}

