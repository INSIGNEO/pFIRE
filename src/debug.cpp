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

#include "debug.hpp"

#include <iomanip>
#include <iostream>

void matrix_dbg_print(const MPI_Comm &comm, const Mat &mat, const std::string &name)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  integer row_lo, row_hi;
  MatGetOwnershipRange(mat, &row_lo, &row_hi);

  MatInfo info;
  MatGetInfo(mat, MAT_LOCAL, &info);

  if (rank == 0)
  {
    std::cout << "Matrix debug info \"" << name << "\":\n"
              << "\tRow range\tmemory\tn_alloc\tn_used\n"
              << std::flush;
  }

  for (int ridx = 0; ridx < size; ridx++)
  {
    if (ridx == rank)
    {
      std::cout << std::setfill(' ') << std::setprecision(2) << "Rank " << rank << ": " << row_lo
                << " - " << row_hi << "\t" << info.memory << "\t" << info.nz_allocated << "\t"
                << info.nz_used << "\n"
                << std::flush;
    }
    MPI_Barrier(comm);
  }
}
