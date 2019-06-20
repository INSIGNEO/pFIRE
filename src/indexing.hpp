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

#ifndef INDEXING_HPP
#define INDEXING_HPP

#include <utility>
#include <vector>

#include "types.hpp"
#include "exceptions.hpp"

/*! Flatten a 3D array index into 1D array location
 *
 */
template <typename inttype>
coord<typename std::enable_if<
    std::is_integral<inttype>::value && !std::is_same<inttype, bool>::value, inttype>::type>
unravel(const inttype& idx, const coord<inttype>& shape)
{
  // This routine only makes sense to use for inttype types
  static_assert(std::is_integral<inttype>::value, "Integral element type required");

  coord<inttype> loc;
  inttype midx = idx;
  for (size_t i = 0; i < shape.size(); i++)
  {
    loc[i] = midx % shape[i];
    midx /= shape[i];
  }

  return loc;
}

/*! Unflatten a 1D array location into 3D array index
 *
 */
template <typename inttype>
typename std::enable_if<std::is_integral<inttype>::value && !std::is_same<inttype, bool>::value,
                        inttype>::type
ravel(const coord<inttype>& loc, const coord<inttype>& shape)
{
  // This routine only makes sense to use for inttype types
  static_assert(std::is_integral<inttype>::value, "Integral element type required");

  size_t end = loc.size() - 1;
  inttype idx = loc[end];
  for (inttype i = end - 1; i >= 0; i--)
  {
    idx *= shape[i];
    idx += loc[i];
  }

  return idx;
}

/*! Calculate offsets for neighbour indices to a given mesh cell using a Von Neumann stencil
 *
 */
template <typename inttype, typename inttype2>
std::vector<typename std::make_signed<typename std::enable_if<
    std::is_integral<inttype>::value && !std::is_same<inttype, bool>::value, inttype>::type>::type>
calculate_von_neumann_offsets(const coord<inttype>& mesh_shape, const inttype2& ndof)
{
  size_t ndim = mesh_shape[2] == 1 ? 2 : 3;

  using inttype_s = typename std::make_signed<inttype>::type;
  std::vector<inttype_s> offsets(2 * ndim + 1, 0);
  size_t centre_idx = ndim;

  inttype offset = 1;
  for (size_t idx = 1; idx <= ndim; idx++)
  {
    offsets[centre_idx + idx] = offset * ndof;
    offsets[centre_idx - idx] = -offset * ndof;

    offset *= mesh_shape[idx - 1];
  }

  return offsets;
}

/*! Calculate offsets for neighbour indices to a mesh cell using a Moore stencil
 *
 */
template <typename inttype, typename inttype2>
std::vector<typename std::make_signed<typename std::enable_if<
    std::is_integral<inttype>::value && !std::is_same<inttype, bool>::value, inttype>::type>::type>
calculate_moore_offsets(const coord<inttype>& mesh_shape, inttype2 ndof)
{
  using inttype_s = typename std::make_signed<inttype>::type;
  std::vector<inttype_s> offsets;

  integer zmin(-1), zmax(1);
  if (mesh_shape[2] == 1)
  {
    zmin = 0;
    zmax = 0;
  }

  for (inttype_s zz = zmin; zz <= zmax; zz++)
  {
    for (inttype_s yy = -1; yy <= 1; yy++)
    {
      for (inttype_s xx = -1; xx <= 1; xx++)
      {
        offsets.push_back(ndof * (xx + mesh_shape[0] * (yy + mesh_shape[1] * zz)));
      }
    }
  }

  return offsets;
}

/*
template <typename T>
bool ranges_overlap(coord<T> a_lo, coord<T> a_hi,
                    coord<T> b_lo, coord<T> b_hi)
{
  for(int i=0; i<3; i++)
  {
    if(b_hi[i] < a_lo[i] || a_hi[i] < b_lo[1])
    {
      return false;
    }
  }

  return true;
}
*/

/*! Given two ranges, return a new range that is the overlap of the two input ranges
 *
 */
template <typename T>
std::pair<coord<T>, coord<T>> get_overlap(coord<T> a_lo, coord<T> a_hi, coord<T> b_lo,
                                            coord<T> b_hi)
{
  coord<T> reslo, reshi;
  for (int i = 0; i < 3; i++)
  {
    reslo[i] = std::max(a_lo[i], b_lo[i]);
    reshi[i] = std::min(a_hi[i], b_hi[i]);
  }

  return std::make_pair(reslo, reshi);
}

template <typename T>
T range_get_nelements(coord<T> lo, coord<T> hi)
{
  T nelem = 1;
  for (int i = 0; i < 3; i++)
  {
    nelem *= hi[i] - lo[i];
  }

  return nelem;
}

template <typename T>
void clamp_location_lo(coord<T>& loc, coord<T> const& min)
{
  for (size_t idx = 0; idx < loc.size(); idx++)
  {
    loc[idx] = loc[idx] >= min[idx] ? loc[idx] : min[idx];
  }
}

template <typename T>
void clamp_location_hi(coord<T>& loc, coord<T> const& max)
{
  for (size_t idx = 0; idx < loc.size(); idx++)
  {
    loc[idx] = loc[idx] <= max[idx] ? loc[idx] : max[idx];
  }
}

template <typename T>
void clamp_location(coord<T>& loc, coord<T> const& min, coord<T> const& max)
{
  for (size_t idx = 0; idx < loc.size(); idx++)
  {
    loc[idx] = loc[idx] >= min[idx] ? loc[idx] : min[idx];
    loc[idx] = loc[idx] <= max[idx] ? loc[idx] : max[idx];
  }
}

template <typename Input1, typename Input2, typename Input3>
floating calculate_scaled_basis_coefficient(Input1 first1, Input1 last1, Input2 first2,
                                            Input3 scale)
{
  floating res = 1;
  auto it1 = first1;
  auto it2 = first2;
  auto scaleit = scale;
  for (; it1 != last1;)
  {
    res *= 1 - std::abs(*it1++ - *it2++) / (*scaleit++);
  }
  if (res > 1.)
  {
    throw InternalError("Coefficient too large");
  }
  return res;
}

template <typename Input1, typename Input2>
floating calculate_unscaled_basis_coefficient(Input1 first1, Input1 last1, Input2 first2)
{
  floating res = 1;
  auto it1 = first1;
  auto it2 = first2;
  for (; it1 != last1;)
  {
    res *= 1 - std::abs(*it1++ - *it2++);
  }
  if (res > 1.)
  {
    throw InternalError("Coefficient too large");
  }
  return res;
}
#endif
