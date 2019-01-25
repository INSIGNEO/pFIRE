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

#ifndef FD_ROUTINES_HPP
#define FD_ROUTINES_HPP

#include <exception>
#include <memory>

#include <petscdmda.h>
#include <petscvec.h>

#include "petsc_helpers.hpp"
#include "types.hpp"

namespace fd
{
Vec_unique gradient_to_global_unique(const DM &dmda, const Vec &localvec, integer dim);

void gradient_existing(const DM &dmda, const Vec &src, Vec &tgt, integer dim);

} // namespace fd

#endif
