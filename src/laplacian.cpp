#include "laplacian.hpp"

#include "indexing.hpp"
#include "petsc_debug.hpp"

Mat_unique build_laplacian_autostride(MPI_Comm comm, intvector shape, integer ndim)
{
  integer matsize = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<integer>());
  matsize *= ndim;
  // get domain info to determine rows
  int rank, num_ranks;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);

  integer rowsize = matsize / num_ranks;
  integer remainder = matsize % num_ranks;
  integer rowstart = rowsize * rank;
  // take care of remainder
  if (rank < remainder)
  {
    rowsize += 1;
    rowstart += rank;
  }
  else
  {
    rowstart += remainder;
  }
  integer rowend = rowstart + rowsize;

  return build_laplacian_matrix(comm, shape, rowstart, rowend, ndim);
}

Mat_unique build_laplacian_matrix(
    MPI_Comm comm, intvector shape, integer startrow, integer endrow, integer ndim)
{
  // total columns == total rows == mask length
  integer n_nodes = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<integer>());
  integer matsize = n_nodes * ndim;

  // Generate laplacian data in CSR format
  intvector idxn, idxm;
  floatvector mdat;
  integer rowptr = 0;
  idxn.push_back(rowptr);
  for (integer gidx = startrow; gidx < endrow; gidx++)
  {
    integer idx = gidx % n_nodes;
    integer ofs = n_nodes * (gidx / n_nodes);
    intvector currloc = unravel(idx, shape);
    integer rowcount = 0;

    for (size_t dim = 0; dim < currloc.size(); dim++)
    {
      currloc[dim] += 1;
      if (currloc[dim] < shape[dim])
      {
        idxm.push_back(ravel(currloc, shape) + ofs);
        mdat.push_back(-0.5);
        rowcount++;
      }
      currloc[dim] -= 2;
      if (currloc[dim] >= 0)
      {
        idxm.push_back(ravel(currloc, shape) + ofs);
        mdat.push_back(-0.5);
        rowcount++;
      }
      currloc[dim] += 1;
    }
    rowptr += rowcount;
    rowptr++;
    idxn.push_back(rowptr);
    idxm.push_back(gidx);
    mdat.push_back(0.5 * rowcount);
  }
  // CSR data consumed directly by PETSc :)
  Mat_unique lapl_mat = create_unique_mat();
  PetscErrorCode perr = MatCreateMPIAIJWithArrays(
      comm, idxn.size() - 1, PETSC_DECIDE, matsize, matsize, idxn.data(), idxm.data(), mdat.data(),
      lapl_mat.get());
  debug_creation(*lapl_mat, "laplacian");
  CHKERRABORT(comm, perr);
  return lapl_mat;
}
