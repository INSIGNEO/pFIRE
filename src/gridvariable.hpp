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

#ifndef GRIDVARIABLE_HPP
#define GRIDVARIABLE_HPP

#include "types.hpp"

#include <algorithm>
#include <numeric>

#include "petscdmda.h"

class GridVariable {
public:
  GridVariable(const intcoord& shape, const integer& ndof, const MPI_Comm& comm);

  //! Returns the MPI communicator of the mesh object
  MPI_Comm comm() const { return _comm; }
  //! Returns the size of the mesh object communicator
  int commsize() const { return _commsize; }
  //! Returns the number of physical dimensions the object has
  integer ndim() const { return _ndim; }
  //! Returns the number of degrees of freedom (vector components) of the grid
  integer ndof() const { return _ndof; }
  //! Returns the shape of the mesh object (3D geometry)
  const intcoord& shape() const { return _shape; }
  //! Returns the offset (bottom corner location) of the locally owned segment of the mesh object
  const intcoord& local_offset() const { return _local_offset; }
  //! Returns the shape of the locally owned segment of the mesh object (3D geometry)
  const intcoord& local_shape() const { return _local_shape; }

  //! Return the min and max values of the grid data
  floatpair minmax() const;

  //! Returns number of nodes owned by the local rank
  integer num_owned_nodes() const { return _local_shape[0] * _local_shape[1] * _local_shape[2]; }
  //! Returns total number of nodes in the mesh (N.B all dofs count as one node)
  integer num_total_nodes() const { return _shape[0] * _shape[1] * _shape[2]; }

  //! Returns coordinates of the lower left vertex of the mesh
  const intcoord& index_min() const { return _index_min; }
  //! Returns coordinates of the upper right vertex of the mesh
  const intcoord& index_max() const { return _index_max; }

  //! Return coordinates of the locally owned area of the grid
  std::pair<intcoord, intcoord> owned_range() const;

  //! Returns a reference to the mesh DMDA
  const DM& dmda() const { return *_dmda; }

  //! Copy vector data in from compatible vector
  void copy_data_from(const GridVariable& src);
  //! Update local ghost cells from global data
  void update_local_vector() const;

  //! Return rank index holding the coordinates of the provided location
  template <typename T>
  typename std::enable_if<(std::is_integral<T>::value || std::is_floating_point<T>::value)
                              && !std::is_same<T, bool>::value,
                          integer>::type
  get_rank_of_loc(coord<T> tgtloc) const;

  std::shared_ptr<GridVariable> copy() const;
  std::shared_ptr<GridVariable> duplicate() const;

  //! Return a duplicate object with the same shape and degrees of freedom - does not copy data
  template <typename T>
  // static std::shared_ptr<typename std::enable_if<std::is_base_of<GridVariable, T>::value,
  // T>::type>
  static std::shared_ptr<T> duplicate(T grid_object);

  //! Return a deep copy of the object including a copy of the mesh data
  template <typename T>
  static std::shared_ptr<typename std::enable_if<std::is_base_of<GridVariable, T>::value, T>::type>
  copy(const T& grid_object);

  //! Scatter all data to rank zero. Use with caution!
  Vec_unique scatter_to_zero() const;

  //! Find MPI rank distribution for a given mesh shape and comm size
  static intcoord find_rank_distribution(const intcoord& datasize, const integer& commsize);

  /*! Returns a reference to the global data vector
   *
   * Protect our data accessors to encourage encapsulation of access to underlying PETSc objects
   * Ideally nobody outside of the gridvariable classes should ever need to know about data
   * layout...
   */
  const Vec& global_vector() const { return *_globalvec; }
  //! Returns a reference to the local data vector
  const Vec& local_vector() const { return *_localvec; }

protected:
  //! Protected copy constructor - see ::copy and ::duplicate for public functions
  GridVariable(const GridVariable& grid_object);

  //! Implicit copies are bad, mm'kay?
  GridVariable() = delete;

  //! Returns a shared pointer to the DMDA
  DM_shared dmda_ptr() const { return _dmda; }

private:
  MPI_Comm _comm;
  int _commsize;

  const intcoord _index_min = {0, 0, 0};
  intcoord _shape;
  intcoord _index_max;
  integer _ndim;
  integer _ndof;
  intcoord _local_shape;
  intcoord _local_offset;

  intcoord _ranks_per_dim;
  intvector2d _rank_edges;

  static const DMBoundaryType _bnd_type = DM_BOUNDARY_GHOSTED;
  static const DMDAStencilType _stencil_type = DMDA_STENCIL_STAR;
  static const integer _stencil_width = 1;

  DM_shared _dmda;
  Vec_unique _globalvec;
  Vec_unique _localvec;

  //! Setup DMDA object
  void initialise_dmda();

  //! Setup data storage vectors
  void initialise_vectors();

  //! Setup additional info variables
  void initialise_additional_info();
};

/*! Namedconstructor companion function to ::copy, copies both metadata and grid data to new
 * object.
 *
 */
template <typename T>
// std::shared_ptr<typename std::enable_if<std::is_base_of<GridVariable, T>::value, T>::type>
std::shared_ptr<T> GridVariable::duplicate(T grid_object)
{
  return std::shared_ptr<T>(new T(grid_object));
}

/*! Namedconstructor companion function to ::duplicate, copies only metadata to new object, not
 * grid data.
 */
template <typename T>
std::shared_ptr<typename std::enable_if<std::is_base_of<GridVariable, T>::value, T>::type>
GridVariable::copy(const T& grid_object)
{
  auto tmp = std::make_shared<T>(grid_object);
  tmp->copy_data_from(grid_object);
  return tmp;
}

/*! Coordinates outside the domain will return the rank holding the nearest valid coordinate
 *
 * \param tgtloc Location of the coordinate
 */
template <typename T>
typename std::enable_if<(std::is_integral<T>::value || std::is_floating_point<T>::value)
                            && !std::is_same<T, bool>::value,
                        integer>::type
GridVariable::get_rank_of_loc(coord<T> tgtloc) const
{
  intcoord rankloc;
  for (int i = 0; i < 3; i++)
  {
    // Upper bound gives first element strictly greater, so subtract one to get lower edge
    auto rankptr = std::upper_bound(_rank_edges[i].cbegin(), _rank_edges[i].cend(), tgtloc[i]);
    // Out of range below
    if (rankptr == _rank_edges[i].cbegin())
    {
      // nudge up to first valid rank
      rankptr++;
    }
    // out of range above
    else if (rankptr == _rank_edges[i].cend())
    {
      // nudge down to last valid rank
      rankptr--;
    }
    rankloc[i] = std::distance(_rank_edges[i].cbegin(), rankptr) - 1;
  }

  // From PETSc code (src/dm/impls/da/da3.c) find ranks are calulated as R = z*NY*NX + y*NX + x for
  // chunk position x,y,z and grid shape NX,NY,NZ
  return (rankloc[2] * _ranks_per_dim[1] + rankloc[1]) * _ranks_per_dim[0] + rankloc[0];
}

#endif // GRIDBARIABLE_HPP
