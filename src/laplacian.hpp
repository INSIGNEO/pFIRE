#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include<numeric>

#include<mpi.h>
#include<petscmat.h>

#include "types.hpp"

Mat_unique build_laplacian_matrix(MPI_Comm comm, intvector shape, integer startrow, integer endrow,
                                  integer ndim);
Mat_unique build_laplacian_autostride(MPI_Comm comm, intvector shape);

#endif

  
