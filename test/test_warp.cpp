#define BOOST_TEST_MODULE workspace
#include "test_common.hpp"

#include<petscdmda.h>

#include "types.hpp"
#include "image.hpp"
#include "map.hpp"
#include "workspace.hpp"

struct envobjs
{
  envobjs() : image(imgshape), map(image, nodespacing), workspace(image, map){}

  intvector imgshape = {5, 9, 7};
  floatvector nodespacing = {3, 3, 3};
  Image image;
  Map map;
  WorkSpace workspace;
};

BOOST_FIXTURE_TEST_SUITE(indexing, envobjs)

  BOOST_AUTO_TEST_CASE(test_warp_image)
  {
    // Set one pixel in image for warping
    floating data = 1.405;
    intvector loc = {3, 4, 5};

    // Access data via dmda and insert value
    {floating ***ptr;
    DMDAVecGetArray(*image.dmda(), *image.global_vec(), &ptr);
    ptr[loc[2]][loc[1]][loc[0]] = data;
    DMDAVecRestoreArray(*image.dmda(), *image.global_vec(), &ptr);
    }

    // Set one dimension of the map to move one pixel in positive x
    
 
  }

BOOST_AUTO_TEST_SUITE_END()
