#include "laplacian.hpp"

//Forward decl for heavy lifting function
Mat dobuild_laplacian_mpi(MPI_Comm comm, intvector shape, integer startrow,
                          integer endrow);


/*Operator build_masked_laplacian(Image Mask){

  integer startrow, endrow, loc_nrows;
  PetscErrorCode ierr;
  ierr = VecGetOwnerShipRange(mask.datavec, &startrow, &endrow);CHKERRRQ(ierr);
  loc_nrows = endrow - startrow;

}*/

Mat create_laplacian_autostride(MPI_Comm comm, intvector shape){
  
  integer n_nodes = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<integer>());

  // get domain info to determine rows
  integer rank, num_ranks, mpi_err;
  mpi_err = MPI_Comm_rank(comm, &rank);
  mpi_err = MPI_Comm_size(comm, &num_ranks);

  integer rowsize = n_nodes / num_ranks;
  integer remainder = n_nodes % num_ranks;
  integer rowstart = rowsize * rank;
  // take care of remainder
  if(rank < remainder){
    rowsize += 1;
    rowstart += rank;
  } else {
    rowstart += remainder;
  }
  integer rowend = rowstart + rowsize;

  Mat laplacian = create_laplacian(comm, shape, rowstart, rowend);

  return laplacian;
}

// Heavy lifting function that actually builds the laplacian, should be called only by one of the
// above user facing functions
Mat create_laplacian(MPI_Comm comm, intvector shape, integer startrow, integer endrow){

  //total columns == total rows == mask length
  integer n_nodes = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<integer>());

  // Generate laplacian data in CSR format
  intvector idxn, idxm;
  floatvector mdat;
  integer rowptr= 0;
  idxn.push_back(rowptr);
  for(integer idx=startrow; idx<endrow; idx++){
    intvector currloc = unravel(idx, shape);
    integer rowcount = 0;

    for(size_t dim=0; dim < currloc.size(); dim++){
      currloc[dim] += 1;
      if(currloc[dim] < shape[dim]){
        idxm.push_back(ravel(currloc, shape));
        mdat.push_back(-0.5);
        rowcount++;
      }
      currloc[dim] -= 2;
      if(currloc[dim] >= 0){
        idxm.push_back(ravel(currloc, shape));
        mdat.push_back(-0.5);
        rowcount++;
      }
      currloc[dim] += 1;
    }
    rowptr += rowcount;
    rowptr++;
    idxn.push_back(rowptr);
    idxm.push_back(idx);
    mdat.push_back(0.5*rowcount);
  }
  // CSR data consumed directly by PETSc :)
  Mat lapl_mat;
  PetscErrorCode perr = MatCreateMPIAIJWithArrays(
    comm, idxn.size()-1, PETSC_DECIDE, n_nodes, n_nodes,
    idxn.data(), idxm.data(), mdat.data(), &lapl_mat);CHKERRABORT(comm, perr);
  return lapl_mat;
  // Wrap it up with a pretty bow and pass it on
//  Operator lapl = Operator(lapl_mat);
//  return lapl;
}
