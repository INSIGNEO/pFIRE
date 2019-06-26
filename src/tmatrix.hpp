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

#ifndef TMATRIX_HPP
#define TMATRIX_HPP

#include <vector>

#include "petscdm.h"
#include "petscvec.h"

#include "types.hpp"

std::pair<Mat_unique, Vec_unique> build_tmat2_and_tmatfm(const Image& fixed, const Image& moved,
                                                         const MapBase& map);

std::pair<floating, floating> calculate_entries(integer row, integer col, integer idof, const MapBase& map,
                                                const std::vector<Vec_unique>& gradient_data,
                                                const Vec& difference_data, MPI_Comm comm);

floatpairvector2d populate_request_vector(intcoordvector2d remote_locs, const DM& dmda,
                                          const Vec& gradient_vec, const Vec& difference_vec);

void communicate_entries_only(integer idof, const MapBase& map, const std::vector<Vec_unique>& gradient_data,
                              const Vec& difference_vec, MPI_Comm comm);

std::vector<Vec_unique> gradients(const Image& fixed, const Image& moved, integer ndof);

void precondition_tmat2(const Mat& tmat2, integer ndof, MPI_Comm _comm);

floating get_condnum_by_poweriter(const Mat& matrix, floating conv_thres, integer max_iter);

floating get_eigenvalue_by_poweriter(const Mat& matrix, floating conv_thres, integer max_iter);

floating diagonal_sum(const Mat& matrix);

void throw_if_idof_inconsistent(MPI_Comm comm, integer idof);

#endif // TMATRIX_HPP
