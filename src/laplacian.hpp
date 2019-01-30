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

#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include <numeric>

#include <mpi.h>
#include <petscmat.h>

#include "types.hpp"

Mat_unique build_laplacian_matrix(
    MPI_Comm comm, intvector shape, integer startrow, integer endrow, integer ndim);
Mat_unique build_laplacian_autostride(MPI_Comm comm, intvector shape);

#endif
