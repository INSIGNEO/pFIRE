#include<numeric>
#include<cstdio>

#include<petscmat.h>
#include<petscvec.h>
#include<petscdmda.h>
#include<petscviewer.h>
#include<mpi.h>

#include "types.hpp"
#include "laplacian.hpp"
#include "image.hpp"
#include "map.hpp"

void mainflow();

int main(int argc, char **argv){
  PetscErrorCode perr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRABORT(PETSC_COMM_WORLD, perr);

  mainflow();

  PetscFinalize();

  return 0;
}

void mainflow(){

  integer nx = 4;
  integer ny = 4;
  integer matsize = 12;

  MPI_Comm &m_comm = PETSC_COMM_WORLD;

  PetscErrorCode perr;

  DM dmda;
  DMDACreate2d(m_comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR,
               nx, ny, PETSC_DECIDE, PETSC_DECIDE, 1, 0, NULL, NULL, &dmda);
  DMSetUp(dmda);
  Vec gvec;
  DMCreateGlobalVector(dmda, &gvec);

  integer x, y, m, n, rank, comm_size;
  MPI_Comm_rank(m_comm, &rank);
  MPI_Comm_size(m_comm, &comm_size);
  DMDAGetCorners(dmda, &x, &y, NULL, &m, &n, NULL);

  PetscSynchronizedPrintf(m_comm, "Rank: %d x: %d - %d y: %d - %d\n", rank, x, x+m, y, y+n);
  PetscSynchronizedFlush(m_comm, PETSC_STDOUT);

  floating **data;
  DMDAVecGetArray(dmda, gvec, &data);
    for(integer ix=x; ix<x+m; ix++)
    {
      for(integer iy=y; iy<y+n; iy++)
      {
        data[iy][ix] = iy*nx + ix; 
      }
    }

  DMDAVecRestoreArray(dmda, gvec, &data);

  Vec tgtvec;
  VecCreateMPI(m_comm, PETSC_DECIDE, nx*ny*2, &tgtvec);
  integer tgtstart, tgtsize;
  VecGetOwnershipRange(gvec, &tgtstart, &tgtsize);
  tgtsize -= tgtstart;
  IS src_is, tgt_is;
  ISCreateStride(m_comm, tgtsize, tgtstart, 1, &src_is); 
  ISCreateStride(m_comm, tgtsize, tgtstart+5, 1, &tgt_is); 
  AO dmdaao;
  DMDAGetAO(dmda, &dmdaao);
  AOApplicationToPetscIS(dmdaao, src_is);

  VecScatter sct;
  VecScatterCreate(gvec, src_is, tgtvec, tgt_is, &sct);
  VecScatterBegin(sct, gvec, tgtvec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(sct, gvec, tgtvec, INSERT_VALUES, SCATTER_FORWARD);

  VecView(tgtvec, PETSC_VIEWER_STDOUT_WORLD);


  Mat bigmat, littlemat;
  MatCreateAIJ(m_comm, PETSC_DECIDE, PETSC_DECIDE, matsize, matsize*2, 2, NULL, 2, NULL, &bigmat);
  MatSetUp(bigmat);
  for(integer idx=0; idx<matsize; idx++)
  {
    MatSetValue(bigmat, idx, idx, (floating)idx, INSERT_VALUES);
    MatSetValue(bigmat, idx, idx+matsize, (floating)idx, INSERT_VALUES);
  }
  MatAssemblyBegin(bigmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(bigmat, MAT_FINAL_ASSEMBLY);

  MatView(bigmat, PETSC_VIEWER_STDOUT_WORLD);

  IS rows, cols;
  ISCreateStride(m_comm, 2, 2*rank, 1, &rows);
  integer cstart, csize;
  MatGetOwnershipRangeColumn(bigmat, &cstart, &csize);
  csize -= cstart;
  ISCreateStride(m_comm, csize, cstart, 1, &cols);

  MatCreateSubMatrix(bigmat, rows, cols, MAT_INITIAL_MATRIX, &littlemat);

  MatView(littlemat, PETSC_VIEWER_STDOUT_WORLD);

}

