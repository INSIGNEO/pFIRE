#define BOOST_TEST_MODULE indexing
#include "test_common.hpp"

#include<petscdmda.h>

#include "types.hpp"
#include "image.hpp"
#include "indexing.hpp"

struct im
{
  im() : image(imgshape)
  {
    BOOST_TEST_MESSAGE("setup image");
    MPI_Comm_size(PETSC_COMM_WORLD, &commsize);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  }

  intvector imgshape = {3, 4, 5};
  Image image;
  int rank;
  int commsize;
};

BOOST_FIXTURE_TEST_SUITE(indexing, im, *butf::enabled())

  BOOST_AUTO_TEST_CASE(test_ravelling)
  {
    integer imlo, imhi;
    VecGetOwnershipRange(*image.global_vec(), &imlo, &imhi);
    // Access data via dmda and insert value
    {integer xlo, ylo, zlo, xw, yw,zw;
    DMDAGetCorners(*image.dmda(), &xlo, &ylo, &zlo, &xw, &yw, &zw);
    floating ***ptr;
    floating const *fptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    VecGetArrayRead(*image.global_vec(), &fptr);
    for(int irank=0; irank<commsize; irank++)
    {
      if(rank == irank){
        for(integer idx=imlo; idx<imhi; idx++)
        {
          intvector loc = unravel(idx, imgshape);
          std::cout << "Idx " << idx << ": " << loc[0] << ", " << loc[1] << ", " << loc[2] << " "
                    << ptr[loc[2]][loc[1]]+loc[0] << " : " << &fptr[idx] << "\n";
          ptr[loc[2]][loc[1]][loc[0]] = idx;
        }
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }
    VecRestoreArrayRead(*image.global_vec(), &fptr);
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);}

    {floating *ptr;
    VecGetArray(*image.global_vec(), &ptr);
    for(int irank=0; irank<commsize; irank++)
    {
      if(rank == irank)
      {
        std::cout << "Rank " << rank << ": ";
        for(integer idx=imlo; idx<imhi; idx++)
        {
          std::cout << ptr[idx-imlo] << " ";
        }
        std::cout << std::endl;
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }
    VecRestoreArray(*image.global_vec(), &ptr);}

  }

  BOOST_AUTO_TEST_CASE(test_unravelling, *butf::disabled())
  {
    integer imlo, imhi;
    VecGetOwnershipRange(*image.global_vec(), &imlo, &imhi);

    // Access data via dmda and insert ravelled value
    integer xlo, ylo, zlo, wx, wy, wz;
    DMDAGetCorners(*image.dmda(), &xlo, &ylo, &zlo, &wx, &wy, &wz);
    floating ***ptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    intvector idxn(3, 0);
    integer& xx=idxn[0];
    integer& yy=idxn[1];
    integer& zz=idxn[2];
    for(xx=xlo; xx<(xlo+wx); xx++)
    {
      for(yy=ylo; yy<(ylo+wy); yy++)
      {
        for(zz=zlo; zz<(zlo+wz); zz++)
        {
          ptr[zz][yy][xx] = ravel(idxn, image.shape());
        }
      }
    }

    {floating ***ptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    for(integer idx=imlo; idx<imhi; idx++)
    {
      intvector loc = unravel(idx, image.shape());
      BOOST_CHECK(ptr[loc[2]][loc[1]][loc[0]] == idx);
    }
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);}
    

  }

BOOST_AUTO_TEST_SUITE_END()
