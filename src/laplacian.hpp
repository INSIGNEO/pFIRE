#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include<numeric>

#include<mpi.h>
#include<petscmat.h>

#include "types.hpp"
#include "indexing.hpp"

Mat create_laplacian(MPI_Comm comm, intvector shape, integer startrow, integer endrow);
Mat create_laplacian_autostride(MPI_Comm comm, intvector shape);

#endif

  
