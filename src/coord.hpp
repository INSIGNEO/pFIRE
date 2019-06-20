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

#ifndef COORD_HPP
#define COORD_HPP

#include <array>
#include <cmath>
#include <ostream>

#include "iterator_routines.hpp"

template <typename classmathtype> 
class coord : public std::array<classmathtype, 3>
{
public: 
  /*! Allow casting between coord types
   *
   * Compiler should unroll loop as this->size() is constexpr
   */
  template <typename coordtype> operator coordtype() const
  {
    coordtype ret;
    for(size_t iidx=0; iidx < this->size(); iidx++)
    {
      ret[iidx] = this->operator[](iidx);
    }
    return ret;
  }

  /*! Coordinate addition
   *
   */
  template <typename mathtype>
  coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())>
  operator+(const coord<mathtype>& a) const
  {
    coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())> ret;
    for(size_t iidx=0; iidx < this->size(); iidx++)
    {
      ret[iidx] = a[iidx] + this->operator[](iidx);
    }
    return ret;
  }

  /*! Coordinate scalar addition
   *
   */
  template <typename mathtype, std::enable_if_t<std::is_arithmetic<mathtype>::value>* = nullptr>
  coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())>
  operator+(mathtype a) const
  {
    coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())> ret;
    for(size_t iidx=0; iidx< this->size(); iidx++)
    {
      ret[iidx] = a + this->operator[](iidx);
    }
    return ret;
  }

  /*! Coordinate subtraction
   *
   */
  template <typename mathtype>
  coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())>
  operator-(const coord<mathtype>& a) const
  {
    coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())> ret;
    for(size_t iidx=0; iidx < this->size(); iidx++)
    {
      ret[iidx] = this->operator[](iidx) - a[iidx] ;
    }
    return ret;
  }

  /*! Coordinate scalar subtraction
   *
   */
  template <typename mathtype, std::enable_if_t<std::is_arithmetic<mathtype>::value>* = nullptr>
  coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())>
  operator-(mathtype a) const
  {
    coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())> ret;
    for(size_t iidx=0; iidx< this->size(); iidx++)
    {
      ret[iidx] = this->operator[](iidx) - a;
    }
    return ret;
  }

  /*! Coordinate scalar multiplication
   *
   */
  template <typename mathtype, std::enable_if_t<std::is_arithmetic<mathtype>::value>* = nullptr>
  coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())>
  operator*(mathtype a) const
  {
    coord<decltype(std::declval<classmathtype>() + std::declval<mathtype>())> ret;
    for(size_t iidx=0; iidx< this->size(); iidx++)
    {
      ret[iidx] = this->operator[](iidx) * a;
    }
    return ret; 
  }

  template <typename mathtype>
  bool operator>(const coord<mathtype>& a) const
  {
    return all_true(this->cbegin(), this->cend(), a.cbegin(), a.cend(), std::greater<>());
  }

  template <typename mathtype>
  bool operator<(const coord<mathtype>& a) const
  {
    return all_true(this->cbegin(), this->cend(), a.cbegin(), a.cend(), std::less<>());
  }

  template <typename mathtype>
  bool operator>=(const coord<mathtype>& a) const
  {
    return all_true(this->cbegin(), this->cend(), a.cbegin(), a.cend(), std::greater_equal<>());
  }

  template <typename mathtype>
  bool operator<=(const coord<mathtype>& a) const
  {
    return all_true(this->cbegin(), this->cend(), a.cbegin(), a.cend(), std::less_equal<>());
  }

  template <typename mathtype>
  bool operator!=(const coord<mathtype>& a) const
  {
    return all_true(this->cbegin(), this->cend(), a.cbegin(), a.cend(), std::not_equal_to<>());
  }

  template <typename mathtype>
  bool operator==(const coord<mathtype>& a) const
  {
    return all_true(this->cbegin(), this->cend(), a.cbegin(), a.cend(), std::equal_to<>());
  }

  coord<classmathtype> reverse() const
  {
    coord<classmathtype> ret;
    ret[0] = this->operator[](2);
    ret[1] = this->operator[](1);
    ret[2] = this->operator[](0);

    return ret;
  }
};

/*! Coordinate scalar addition
 *
 */
template <typename mathtype1, typename mathtype2,
  std::enable_if_t<std::is_arithmetic<mathtype1>::value>* = nullptr>
coord<decltype(std::declval<mathtype1>() + std::declval<mathtype2>())>
operator+(mathtype1 a, const coord<mathtype2>& b)
{
  return b + a;
}

/*! Coordinate scalar subtraction
 *
 */

template <typename mathtype1, typename mathtype2,
  std::enable_if_t<std::is_arithmetic<mathtype1>::value>* = nullptr>
coord<decltype(std::declval<mathtype1>() - std::declval<mathtype2>())>
operator-(mathtype1 a, const coord<mathtype2>& b)
{
  return - (b - a);
}

/*! Coordinate scalar multiplication
 *
 */

template <typename mathtype1, typename mathtype2,
  std::enable_if_t<std::is_arithmetic<mathtype1>::value>* = nullptr>
coord<decltype(std::declval<mathtype1>() + std::declval<mathtype2>())>
operator*(mathtype1 a, const coord<mathtype2>& b)
{
  return b * a;
}

template <typename mathtype>
coord<mathtype> floor(coord<mathtype> a)
{
  coord<mathtype> res;
  for(size_t iidx=0; iidx < a.size(); iidx++)
  {
    res[iidx] = std::floor(a[iidx]);
  }
  return res;
}

/*! Overload << operator for printing
 *
 */
template <typename mathtype>
std::ostream& operator<<(std::ostream& os, coord<mathtype> a)
{
  return os << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
}
#endif //COORD_HPP
