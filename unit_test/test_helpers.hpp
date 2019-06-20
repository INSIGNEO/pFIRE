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

#ifndef TEST_HELPERS_HPP
#define TEST_HELPERS_HPP

#include <type_traits>
#include <cmath>

template <typename ftype>
typename std::enable_if<std::is_floating_point<ftype>::value, bool>::type
fcmp(ftype a, ftype b, ftype rtol, ftype atol)
{
  return std::fabs(a-b) <= atol + rtol*std::fabs(b);
}

#endif //TEST_HELPERS_HPP
