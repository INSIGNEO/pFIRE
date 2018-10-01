#include<numeric>

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

  PetscErrorCode perr;

  intvector shape = {100,100,100};
 
  PetscPrintf(PETSC_COMM_WORLD,"Creating laplacian\n");
  Mat lapl = create_laplacian_autostride(PETSC_COMM_WORLD, shape);

  PetscPrintf(PETSC_COMM_WORLD,"Saving laplacian\n");
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, "lapl.dat", FILE_MODE_WRITE, &viewer);
  MatView(lapl, viewer);
  PetscViewerDestroy(&viewer);

  PetscPrintf(PETSC_COMM_WORLD, "Creating testimg\n");
  intvector imshape = {50, 50, 1};
  Image testimg(imshape);

  PetscPrintf(PETSC_COMM_WORLD, "Calculating gradients\n");
  auto gradx = testimg.gradient(0);
  auto grady = testimg.gradient(1);

  DMView(*(testimg.dmda()), PETSC_VIEWER_STDOUT_WORLD);
//  Image gradz = testimg.gradient(2);

  PetscPrintf(PETSC_COMM_WORLD, "Creating test matrices");
  Mat testmat;
  MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 5, 5, 3, NULL, 3, NULL, &testmat);
  MatSetUp(testmat);
  MatSetValue(testmat, 0, 0, 1, INSERT_VALUES);
  MatSetValue(testmat, 1, 0, 1, INSERT_VALUES);
  MatSetValue(testmat, 1, 1, 1, INSERT_VALUES);
  MatSetValue(testmat, 2, 0, 1, INSERT_VALUES);
  MatSetValue(testmat, 2, 2, 1, INSERT_VALUES);
  MatSetValue(testmat, 4, 3, 1, INSERT_VALUES);
  MatAssemblyBegin(testmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(testmat, MAT_FINAL_ASSEMBLY);
  MatView(testmat, PETSC_VIEWER_STDOUT_WORLD);

  Vec testvec;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 5, &testvec);
  VecSetValue(testvec, 0, 1, INSERT_VALUES);
  VecSetValue(testvec, 1, 2, INSERT_VALUES);
  VecSetValue(testvec, 2, 3, INSERT_VALUES);
  VecSetValue(testvec, 3, 4, INSERT_VALUES);
  VecSetValue(testvec, 4, 5, INSERT_VALUES);
  VecAssemblyBegin(testvec);
  VecAssemblyEnd(testvec);

  MatDiagonalScale(testmat, testvec, NULL);
  MatView(testmat, PETSC_VIEWER_STDOUT_WORLD);

  Mat nest;
  Mat submats[2] = {testmat, testmat};
  MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 1, NULL, submats, &nest);
  MatAssemblyBegin(nest, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(nest, MAT_FINAL_ASSEMBLY);
  MatView(nest, PETSC_VIEWER_STDOUT_WORLD);

  Vec nesttestvec, nesttestoutvec;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 10, &nesttestvec);
  VecDuplicate(nesttestvec, &nesttestoutvec);
  VecSet(nesttestvec, 1.0);
  Mat res;
  perr = MatMatTransposeMult(nest, nest, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &res);CHKERRABORT(PETSC_COMM_WORLD, perr);

  MatView(res, PETSC_VIEWER_STDOUT_WORLD);


}

