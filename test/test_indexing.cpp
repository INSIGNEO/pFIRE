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
    floating data = 1.023;

    // Access data via dmda and insert value
    {floating ***ptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    ptr[3][2][1] = data;
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);
    }
    
    // Access data via vector and retrieve from same place.
    intvector loc = {1,2,3};
    int idx = ravel(loc, imgshape); 
    std::printf("idx: %i\n", idx);
    floating *ptr;
    VecGetArray(*image.global_vec(), &ptr);
    BOOST_CHECK(ptr[idx] == data);
    VecRestoreArray(*image.global_vec(), &ptr);

  }

  BOOST_AUTO_TEST_CASE(test_unravelling)
  {
    floating data = 1.405;
    integer idx = 231;

    // Access data via vector and retrieve from same place.
    std::printf("idx: %i\n", idx);
    floating *ptr;
    VecGetArray(*image.global_vec(), &ptr);
    ptr[idx] = data;
    VecRestoreArray(*image.global_vec(), &ptr);

    // Access data via dmda and insert value
    {floating ***ptr;
    intvector loc = unravel(idx, imgshape);
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    BOOST_CHECK(ptr[loc[2]][loc[1]][loc[0]] == data);
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);
    }
    

  }


BOOST_AUTO_TEST_SUITE_END()
