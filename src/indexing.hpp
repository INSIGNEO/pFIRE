#ifndef INDEXING_HPP
#define INDEXING_HPP

#include <vector>

#include "types.hpp"

template <typename inttype>
std::vector<inttype> unravel(const inttype& idx, const std::vector<inttype>& shape)
{
  // This routine only makes sense to use for inttype types
  static_assert(std::is_integral<inttype>::value, "Integral element type required");

  std::vector<inttype> loc(shape.size());
  inttype midx = idx;
  for (inttype i = 0; i < (inttype)shape.size(); i++)
  {
    loc[i] = midx % shape[i];
    midx /= shape[i];
  }

  return loc;
}

template <typename inttype>
inline inttype ravel(const std::vector<inttype>& loc, const std::vector<inttype>& shape)
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

template <typename inttype>
inline inttype idx_cmaj_to_rmaj(const inttype& idx, const std::vector<inttype>& shape)
{
  inttype z = idx / (shape[0] * shape[1]);
  inttype y = (idx / shape[0]) % shape[1];
  inttype x = idx % shape[0];

  return x * shape[1] * shape[2] + y * shape[2] + z;
}


#endif
