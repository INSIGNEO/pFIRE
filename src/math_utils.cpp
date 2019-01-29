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
