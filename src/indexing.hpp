#ifndef INDEXING_HPP
#define INDEXING_HPP

#include<vector>

#include "types.hpp"

inline intvector unravel(const integer& idx, const intvector& shape)
{
  intvector loc(shape.size());
  integer midx = idx;
  for(integer i=0; i < (integer)shape.size(); i++){
    loc[i] = midx%shape[i];
    midx /= shape[i];
  }
  
  return loc;
}

inline integer ravel(const intvector& loc, const intvector& shape)
{
  size_t end = loc.size() - 1;
  integer idx = loc[end];
  for(integer i=end-1; i >= 0; i--){
    idx *= shape[i];
    idx += loc[i];
  }

  return idx;
}

#endif
