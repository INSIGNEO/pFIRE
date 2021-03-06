#define BOOST_TEST_MODULE indexing
#include "test_common.hpp"

#include<petscdmda.h>
#include<petscerror.h>

#include "types.hpp"
#include "image.hpp"
#include "indexing.hpp"
#include "fd_routines.hpp"

struct im
{
  im() : image(imgshape) {BOOST_TEST_MESSAGE("setup image");}

  intvector imgshape = {3, 3, 3};
  Image image;
};

BOOST_FIXTURE_TEST_SUITE(indexing, im)

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

    Vec_unique grad = image.gradient(0);

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
          BOOST_REQUIRE(ptr[zz][yy][xx] == 1);
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

    Vec_unique grad = image.gradient(1);

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
          BOOST_REQUIRE(ptr[zz][yy][xx] == 1);
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

    Vec_unique grad = image.gradient(2);

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
          BOOST_REQUIRE(ptr[zz][yy][xx] == 3);
        }
      }
    }
    perr = DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);CHKERRXX(perr);
    }//gradient checking enclosure
  }

BOOST_AUTO_TEST_SUITE_END()
