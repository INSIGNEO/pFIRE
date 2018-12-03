#include "debug.hpp"

#include<iomanip>
#include<iostream>

void matrix_dbg_print(const MPI_Comm &comm, const Mat &mat, const std::string &name)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  integer row_lo, row_hi;
  MatGetOwnershipRange(mat, &row_lo, &row_hi);

  MatInfo info;
  MatGetInfo(mat, MAT_LOCAL, &info);

  
  if(rank == 0)
  {
    std::cout << "Matrix debug info \"" << name << "\":\n" 
              << "Row range\tmemory\tnz_allocated\tnz_used"
              << std::flush;
  }

  for(int ridx = 0; ridx < size; ridx++)
  {
    if(ridx == rank)
    {
      std::cout << std::setfill(' ')
                << "Rank " << rank << ": " << row_lo << " - " << row_hi << "\t "
                << info.memory << "\t" << info.nz_allocated << "\t"
                << info.nz_used << "\n" << std::flush;
    }
    MPI_Barrier(comm);
  }
}
