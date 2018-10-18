#define BOOST_TEST_MODULE indexing
#include "test_common.hpp"

#include<petscdmda.h>

#include "types.hpp"
#include "image.hpp"
#include "indexing.hpp"

struct im
{
  im() : image(imgshape) {BOOST_TEST_MESSAGE("setup image");}

  intvector imgshape = {10, 13,11};
  Image image;
};

BOOST_FIXTURE_TEST_SUITE(indexing, im)

  BOOST_AUTO_TEST_CASE(test_ravelling)
  {
    integer imlo, imhi;
    VecGetOwnershipRange(*image.global_vec(), &imlo, &imhi);
    // Access data via dmda and insert value
    {floating ***ptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    for(integer idx=imlo; idx<imhi; idx++)
    {
      intvector loc = unravel(idx, imgshape);
      ptr[loc[2]][loc[1]][loc[0]] = idx;
    }
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);}

    {floating *ptr;
    VecGetArray(*image.global_vec(), &ptr);
    for(integer idx=imlo; idx<imhi; idx++)
    {
      BOOST_CHECK(ptr[idx-imlo] == idx);
    }
    VecRestoreArray(*image.global_vec(), &ptr);}

  }

  BOOST_AUTO_TEST_CASE(test_unravelling)
  {
    integer imlo, imhi;
    VecGetOwnershipRange(*image.global_vec(), &imlo, &imhi);
    // Access data via dmda and insert value

    {floating *ptr;
    VecGetArray(*image.global_vec(), &ptr);
    for(integer idx=imlo; idx<imhi; idx++)
    {
      ptr[idx-imlo] = idx;
    }
    VecRestoreArray(*image.global_vec(), &ptr);}

    {floating ***ptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    for(integer idx=imlo; idx<imhi; idx++)
    {
      intvector loc = unravel(idx, imgshape);
      BOOST_CHECK(ptr[loc[2]][loc[1]][loc[0]] == idx);
    }
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);}
    

  }

BOOST_AUTO_TEST_SUITE_END()
