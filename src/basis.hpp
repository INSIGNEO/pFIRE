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

#ifndef BASIS_HPP
#define BASIS_HPP

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <stdexcept>

#include "types.hpp"

Mat_unique build_basis_matrix(
    MPI_Comm comm, const intvector& src_shape, const intvector& tgt_shape,
    const floatvector& scalings, const floatvector& offsets, uinteger ndim, uinteger tile_dim);

Mat_unique build_warp_matrix(
    MPI_Comm comm, const intvector& img_shape, uinteger ndim,
    const std::vector<Vec*>& displacements);

template <
    class Input1, class Input2, class Rtype = typename std::iterator_traits<Input1>::value_type>
Rtype calculate_basis_coefficient(Input1 first1, Input1 last1, Input2 first2)
{
  Rtype res = 1.;
  auto it1 = first1;
  auto it2 = first2;
  for (; it1 != last1;)
  {
    floating diff = std::abs(*it1++ - *it2++);
    if(diff >= 1.)
    {
      return 0.;
    }
    res *= 1. - diff;
  }
  return res;
}

inline floating clamp_to_edge(floating idx, integer dimsize)
{
  return (idx < 0.) ? 0. : ((idx > dimsize - 2) ? dimsize - 1 : idx);
}

#endif
