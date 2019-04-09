#define BOOST_TEST_MODULE image
#include "test_common.hpp"

#include<petscdmda.h>
#include<petscerror.h>

#include "types.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "fd_routines.hpp"

struct im
{
  im() : image(imgshape)
  {
    MPI_Comm_size(PETSC_COMM_WORLD, &commsize);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  }

  int rank;
  int commsize;

  intvector imgshape = {30, 30, 30};
  Image image;
};

BOOST_FIXTURE_TEST_SUITE(image, im)

  BOOST_AUTO_TEST_CASE(test_xgradient)
  {
    PetscErrorCode perr;

    // Access data via dmda and insert value
    integer xlo, xhi, ylo, yhi, zlo, zhi;
    perr = DMDAGetCorners(*image.dmda(), &xlo, &ylo, &zlo, &xhi, &yhi, &zhi);CHKERRXX(perr);
    xhi += xlo;
    yhi += ylo;
    zhi += zlo;
    {//image setting enclosure
    floating ***ptr;
    perr = DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    for(integer xx=xlo; xx<xhi; xx++)
    {
      for(integer yy=ylo; yy<yhi; yy++)
      {
        for(integer zz=zlo; zz<zhi; zz++)
        {
          ptr[zz][yy][xx] = xx;
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//image setting enclosure

    image.update_local_from_global();
    Vec_unique grad = fd::gradient_to_global_unique(*image.dmda(), *image.local_vec(), 0);

    {//gradient checking enclosure
    floating ***ptr;
    perr = DMDAVecGetArray(*image.dmda(), *grad, &ptr);CHKERRXX(perr);
    xlo = (xlo >= imgshape[0]-1) ? imgshape[0]-1 : ((xlo < 1) ? 1 : xlo); // clamp x range to between 1 and n-1
    xhi = (xhi <= 1) ? 1 : ((xhi > imgshape[0]-1) ? imgshape[0]-1 : xhi); // but ensure at least 1 row checked
    for(integer xx=xlo; xx<xhi; xx++)
    {
      for(integer yy=ylo; yy<yhi; yy++)
      {
        for(integer zz=zlo; zz<zhi; zz++)
        {
          BOOST_CHECK(ptr[zz][yy][xx] == 1);
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//gradient checking enclosure
  }


  BOOST_AUTO_TEST_CASE(test_ygradient)
  {
    PetscErrorCode perr;

    // Access data via dmda and insert value
    integer xlo, xhi, ylo, yhi, zlo, zhi;
    perr = DMDAGetCorners(*image.dmda(), &xlo, &ylo, &zlo, &xhi, &yhi, &zhi);CHKERRXX(perr);
    xhi += xlo;
    yhi += ylo;
    zhi += zlo;
    {//image setting enclosure
    floating ***ptr;
    perr = DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    for(integer xx=xlo; xx<xhi; xx++)
    {
      for(integer yy=ylo; yy<yhi; yy++)
      {
        for(integer zz=zlo; zz<zhi; zz++)
        {
          ptr[zz][yy][xx] = yy;
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//image setting enclosure

    image.update_local_from_global();
    Vec_unique grad = fd::gradient_to_global_unique(*image.dmda(), *image.local_vec(), 1);

    {//gradient checking enclosure
    floating ***ptr;
    perr = DMDAVecGetArray(*image.dmda(), *grad, &ptr);CHKERRXX(perr);
    ylo = (ylo >= imgshape[1]-1) ? imgshape[1]-1 : ((ylo < 1) ? 1 : ylo); // clamp x range to between 1 and n-1
    yhi = (yhi <= 1) ? 1 : ((yhi > imgshape[1]-1) ? imgshape[1]-1 : yhi); // but ensure at least 1 row checked
    for(integer xx=xlo; xx<xhi; xx++)
    {
      for(integer yy=ylo; yy<yhi; yy++)
      {
        for(integer zz=zlo; zz<zhi; zz++)
        {
          BOOST_CHECK(ptr[zz][yy][xx] == 1);
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//gradient checking enclosure
  }


  BOOST_AUTO_TEST_CASE(test_zgradient)
  {
    PetscErrorCode perr;

    // Access data via dmda and insert value
    integer xlo, xhi, ylo, yhi, zlo, zhi;
    perr = DMDAGetCorners(*image.dmda(), &xlo, &ylo, &zlo, &xhi, &yhi, &zhi);CHKERRXX(perr);
    xhi += xlo;
    yhi += ylo;
    zhi += zlo;
    {//image setting enclosure
    floating ***ptr;
    perr = DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    for(integer xx=xlo; xx<xhi; xx++)
    {
      for(integer yy=ylo; yy<yhi; yy++)
      {
        for(integer zz=zlo; zz<zhi; zz++)
        {
          ptr[zz][yy][xx] = 3*zz;
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//image setting enclosure

    image.update_local_from_global();
    Vec_unique grad = fd::gradient_to_global_unique(*image.dmda(), *image.local_vec(), 2);

    {//gradient checking enclosure
    floating ***ptr;
    perr = DMDAVecGetArray(*image.dmda(), *grad, &ptr);CHKERRXX(perr);
    zlo = (zlo >= imgshape[2]-1) ? imgshape[2]-1 : ((zlo < 1) ? 1 : zlo); // clamp x range to between 1 and n-1
    zhi = (zhi <= 1) ? 1 : ((zhi > imgshape[2]-1) ? imgshape[2]-1 : zhi); // but ensure at least 1 row checked
    for(integer xx=xlo; xx<xhi; xx++)
    {
      for(integer yy=ylo; yy<yhi; yy++)
      {
        for(integer zz=zlo; zz<zhi; zz++)
        {
          BOOST_CHECK(ptr[zz][yy][xx] == 3);
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//gradient checking enclosure
  }

  BOOST_AUTO_TEST_CASE(test_get_rank_of_loc)
  {
    // Access data via dmda and insert value
    integer xlo, xhi, ylo, yhi, zlo, zhi;
    PetscErrorCode perr;
    perr = DMDAGetCorners(*image.dmda(), &xlo, &ylo, &zlo, &xhi, &yhi, &zhi);CHKERRXX(perr);
    xhi += xlo;
    yhi += ylo;
    zhi += zlo;

    floatvector xlocs = {double(xlo), xhi-1.1};
    floatvector ylocs = {double(ylo), yhi-1.1};
    floatvector zlocs = {double(zlo), zhi-1.1};

    for(const auto &xx: xlocs)
    {
      for(const auto &yy: ylocs)
      {
        for(const auto &zz: zlocs)
        {
          std::cout << rank;
          integer rl = image.get_rank_of_loc({xx, yy, zz});
          BOOST_CHECK(rl == rank);
        }
      }
    }

  }

BOOST_AUTO_TEST_SUITE_END()
