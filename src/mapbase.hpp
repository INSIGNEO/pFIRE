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

#ifndef MAPBASE_HPP
#define MAPBASE_HPP

#include "gridvariable.hpp"
#include "mask.hpp"
#include "types.hpp"

/*! "Serial" Map Class
 *
 * These classes are a bit hackish due to the fact that PETSc does not support DMDAs being defined on more
 * processors than there are nodes in the grid.  Therefore, we have both a Parallel and a "Serial" class.
 * For the serial clas, while the DMDA is defined on MPI_COMM_SELF, the actual map data is
 * still distributed to enable efficient parallel computation of T^tT and T^t(f-m). 
 *
 *  Because of this, we define some extra functions to allow the solver to interact with a MapBase class
 *  while being oblivious to any nasty implementation details lurking under the hood...
 *
 *  Unfortunately this makes interacting with a MapBase object through the inherited global_vector() and
 *  local_vector() methods a bit risky as depending on the true nature of this object, the vector will either
 *  be defined on COMM_WORLD or COMM_SELF and so mathematical sanity may rapidly vanish...
 *
 * */
class MapBase: public GridVariable {
public:
  MapBase(const intcoord& shape, const integer& ndim, const intcoord& node_spacing, const Mask& mask,
          const MPI_Comm& comm);

  const floatvector2d& node_locs() const { return _node_locations; };
  const intcoord& spacing() const { return _node_spacing; }
  const Mask& mask() const { return _mask; };
  const integer& stencilpoints() const { return _stencilpoints; };

  std::pair<intcoord, intcoord> get_pixel_neighbourhood(const intcoord& map_node) const;

  floatcoord coord_from_index(const intcoord& index) const;

  //! Return the node at the lower corner of the cell containing loc
  template<typename mathtype>
  intcoord get_cell_corner(coord<mathtype> loc) const;

  //! Given a mask and node spacing, create a map of appropriate size and shape.
  static std::unique_ptr<MapBase> make_map_for_mask(const Mask& mask,
                                                    const intcoord& node_spacing);

  static intcoord calculate_map_shape(const intcoord& mask_shape, const intcoord& target_spacing);

  floatcoord lower_corner() const {return this->coord_from_index(this->index_min());}
  floatcoord upper_corner() const {return this->coord_from_index(this->index_max());}

  virtual bool parallel_layout() const=0;

  virtual const Vec& displacement_vector() const=0;

  const Mat& laplacian2() const;

  virtual void interpolate_from(const MapBase& previous_map)=0;

  virtual void update(const Vec& delta)=0;

  floatcoordvector map_local_coordinates(const intcoordvector& loc) const;

  static const intcoordvector node_offset_list_2d;
  static const intcoordvector node_offset_list_3d;

  void map_print() const;

protected:

private:
  const Mask& _mask;
  const integer _stencilpoints = 3;
  intcoord _node_spacing;
  floatvector2d _node_locations;
  floatvector _offsets;

  mutable Mat_unique _laplacian2;

  void calculate_node_locations();

  void build_laplacian2() const;

  static MPI_Comm validate_map_comm(const MPI_Comm& comm, const Mask& mask);
};

//! Return the node at the lower corner of the cell containing loc
template<typename mathtype>
intcoord MapBase::get_cell_corner(coord<mathtype> loc) const
{
  intcoord node_loc{0, 0, 0};
  for (integer idim = 0; idim < this->ndim(); idim++)
  {
    auto iptr = std::upper_bound(_node_locations[idim].cbegin(), _node_locations[idim].cend(), loc[idim]);
    if (iptr == _node_locations[idim].cend())
    {
      node_loc[idim] = _node_locations[idim].size() - 2;
      continue;
    }
    if (iptr == _node_locations[idim].cbegin())
    {
      node_loc[idim] = 0;
      continue;
    }
    node_loc[idim] = std::distance(_node_locations[idim].cbegin(), iptr) - 1;
  }

  return node_loc;
}

#endif // MAPBASE_HPP
