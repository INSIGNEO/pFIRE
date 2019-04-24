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

#include "math_utils.hpp"

#include <petscdmda.h>

void quadratic_from_points(floating x_1, floating x_2, floating x_3,
                           floating y_1, floating y_2, floating y_3,
                           floating &a, floating &b, floating &c)
{

  a = y_1/((x_1-x_2)*(x_1-x_3)) + y_2/((x_2-x_1)*(x_2-x_3)) + y_3/((x_3-x_1)*(x_3-x_2));

  b = -y_1*(x_2+x_3)/((x_1-x_2)*(x_1-x_3))
      -y_2*(x_1+x_3)/((x_2-x_1)*(x_2-x_3))
      -y_3*(x_1+x_2)/((x_3-x_1)*(x_3-x_2));

  c = y_1*x_2*x_3/((x_1-x_2)*(x_1-x_3))
    + y_2*x_1*x_3/((x_2-x_1)*(x_2-x_3))
    + y_3*x_1*x_2/((x_3-x_1)*(x_3-x_2));
}

void quadratic_vertex(floating a, floating b, floating c, floating &x, floating &y)
{
  x = -b/(2*a);

  y = a*x*x + b*x + c;
}

void binarize_vector(Vec &vec, floating thres)
{
  // Binarize image using threshold in fractional range 0-1 of image range.
  integer localsize;
  PetscErrorCode perr = VecGetLocalSize(vec, &localsize);
  CHKERRXX(perr);

  floating vecmax;
  perr = VecMax(vec, nullptr, &vecmax);
  CHKERRXX(perr);
  perr = VecScale(vec, 1/vecmax);
  CHKERRXX(perr);

  floating *data;
  perr = VecGetArray(vec, &data);
  CHKERRXX(perr);

  for(integer idx = 0; idx < localsize; idx++)
  {
    data[idx] = data[idx] > thres ? 1.0 : 0.0;
  }

  perr = VecRestoreArray(vec, &data);
  CHKERRXX(perr);
}


void dilate_dmda(DM &dmda, Vec &vec)
{
  // Acquire locally owned segment
  integer x_lo, y_lo, z_lo, x_hi, y_hi, z_hi;
  PetscErrorCode perr = DMDAGetCorners(dmda, &x_lo, &y_lo, &z_lo, &x_hi, &y_hi, &z_hi);
  CHKERRXX(perr);
  x_hi += x_lo;
  y_hi += y_lo;
  z_hi += z_lo;

  Vec_unique dilation_vec = create_unique_vec();
  perr = DMCreateLocalVector(dmda, dilation_vec.get());
  CHKERRXX(perr);
  perr = VecSet(*dilation_vec, 0);

  floating ***vecdata, ***dildata;
  perr = DMDAVecGetArray(dmda, vec, &vecdata);
  CHKERRXX(perr);
  perr = DMDAVecGetArray(dmda, *dilation_vec, &dildata);
  CHKERRXX(perr);

  for(integer x_idx = x_lo; x_idx < x_hi; x_idx++)
  {
    for(integer y_idx = y_lo; y_idx < y_hi; y_idx++)
    {
      for(integer z_idx = z_lo; z_idx < z_hi; z_idx++)
      {
        if(vecdata[z_idx][y_idx][x_idx] > 0)
        {
          dildata[z_idx+1][y_idx][x_idx] += vecdata[z_idx][y_idx][x_idx];
          dildata[z_idx-1][y_idx][x_idx] += vecdata[z_idx][y_idx][x_idx];
          dildata[z_idx][y_idx+1][x_idx] += vecdata[z_idx][y_idx][x_idx];
          dildata[z_idx][y_idx-1][x_idx] += vecdata[z_idx][y_idx][x_idx];
          dildata[z_idx][y_idx][x_idx+1] += vecdata[z_idx][y_idx][x_idx];
          dildata[z_idx][y_idx][x_idx-1] += vecdata[z_idx][y_idx][x_idx];
        }
      }
    }
  }

  perr = DMDAVecRestoreArray(dmda, vec, &vecdata);
  CHKERRXX(perr);
  perr = DMDAVecRestoreArray(dmda, *dilation_vec, &dildata);
  CHKERRXX(perr);

  perr = DMLocalToGlobalBegin(dmda, *dilation_vec, ADD_VALUES, vec);
  CHKERRXX(perr);
  perr = DMLocalToGlobalEnd(dmda, *dilation_vec, ADD_VALUES, vec);
  CHKERRXX(perr);
  
}
