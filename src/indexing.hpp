#ifndef INDEXING_HPP
#define INDEXING_HPP

#include<vector>

#include "types.hpp"

inline intvector unravel(const integer& idx, const intvector& shape)
{
  intvector loc(shape.size());
  integer midx = idx;
  for(integer i=shape.size()-1; i >= 0; i--){
    loc[i] = midx%shape[i];
    midx /= shape[i];
  }
  
  return loc;
}

inline integer ravel(const intvector& loc, const intvector& shape)
{
  integer idx = loc[0];
  for(integer i=1; i < (integer)shape.size(); i++){
    idx *= shape[i];
    idx += loc[i];
  }

  return idx;
}


#endif
