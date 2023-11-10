//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//

#define BOOST_TEST_MODULE workspace
#include "test_common.hpp"

#include <petscdmda.h>

#include "image.hpp"
#include "mapbase.hpp"
#include "mask.hpp"
#include "types.hpp"
#include "test_helpers.hpp"
#include "xdmfwriter.hpp"

class envobjs {
public:
  envobjs()
    : image(std::make_unique<Image>(imgshape, MPI_COMM_WORLD)), mask(Mask::create_filled_mask(*image)),
      map(MapBase::make_map_for_mask(*mask, nodespacing))
  {
  }

  const floating testdata = 1.0;
  const integer offset = 5;
  const intcoord imgshape = {50, 50, 50};
  const intcoord nodespacing = {30, 30, 30};
  const intcoord block_shape = {30,30,30};
  const intcoord block_offset = {10,10,10};
  std::shared_ptr<Image> image;
  std::shared_ptr<Mask> mask;
  std::shared_ptr<MapBase> map;
};

BOOST_FIXTURE_TEST_SUITE(indexing, envobjs)

BOOST_AUTO_TEST_CASE(uniform_warp)
{

  for(integer curr_dof = 0; curr_dof < map->ndof(); curr_dof++)
  {
    std::cout << "Iteration " << curr_dof << "\n";
    // reset map
    map = MapBase::make_map_for_mask(*mask, nodespacing);
    image = std::make_unique<Image>(imgshape, MPI_COMM_WORLD);

    intcoord image_local_lo = image->local_offset();
    intcoord image_local_hi = image_local_lo + image->local_shape();
    // Access data via dmda and insert value
    floating*** image_data;
    DMDAVecGetArray(image->dmda(), image->data_vector(), &image_data);
    intcoord curr_loc;
    auto& xx = curr_loc[0];
    auto& yy = curr_loc[1];
    auto& zz = curr_loc[2];
    for(zz=image_local_lo[2]; zz < image_local_hi[2]; zz++)
    {
      for(yy=image_local_lo[1]; yy < image_local_hi[1]; yy++)
      {
        for(xx=image_local_lo[0]; xx < image_local_hi[0]; xx++)
        {
          if(curr_loc >= block_offset && curr_loc < (block_offset + block_shape))
          {
            image_data[zz][yy][xx] = testdata;
          }
        }
      }
    }
    DMDAVecRestoreArray(image->dmda(), image->data_vector(), &image_data);
    image->update_local_vector();

    Vec_unique delta = create_unique_vec();
    VecDuplicate(map->displacement_vector(), delta.get());
    // Set one dimension of the map to move one pixel in positive direction
    VecSet(*delta, 0.);
    floating* delta_data;
    integer veclo, vechi;
    VecGetOwnershipRange(*delta, &veclo, &vechi);
    VecGetArray(*delta, &delta_data);
    for(integer iidx=veclo; iidx<vechi; iidx++)
    {
      if(iidx%map->ndof() == curr_dof)
      {
        delta_data[iidx-veclo] = offset;
      }
    }
    VecRestoreArray(*delta, &delta_data);

    map->update(*delta);
    map->update_local_vector();

    floating**** map_data;
    DMDAVecGetArrayDOF(map->dmda(), map->global_vector(), &map_data);
    intcoord map_local_lo = map->local_offset();
    intcoord map_local_hi = map_local_lo + map->local_shape();
    for(integer idof=0; idof<4; idof++)
      {
        floating cmp = (idof == curr_dof) ? offset : 0.;
        for(integer zz=map_local_lo[2]; zz<map_local_hi[2]; zz++)
        {
          for(integer yy=map_local_lo[1]; yy<map_local_hi[1]; yy++)
          {
            for(integer xx=map_local_lo[0]; xx<map_local_hi[0]; xx++)
            {
              BOOST_REQUIRE(fcmp((double) map_data[zz][yy][xx][idof],(double) cmp, 1e-10, 1e-10));
            }
          }
        }
      }
    DMDAVecRestoreArrayDOF(map->dmda(), map->global_vector(), &map_data);

    auto warped = image->warp(*map);

    DMDAVecGetArray(warped->dmda(), warped->data_vector(), &image_data);

    intcoord warp_offset({0, 0, 0});
    if(curr_dof > 0)
    {
      warp_offset[curr_dof-1] = offset;
    }

    DMDAVecGetArray(warped->dmda(), warped->data_vector(), &image_data);
    for(zz=image_local_lo[2]; zz < image_local_hi[2]; zz++)
    {
      for(yy=image_local_lo[1]; yy < image_local_hi[1]; yy++)
      {
        for(xx=image_local_lo[0]; xx < image_local_hi[0]; xx++)
        {
          if((curr_loc-warp_offset) >= block_offset && (curr_loc-warp_offset) < (block_offset + block_shape))
          {
            BOOST_REQUIRE(fcmp((double) image_data[zz][yy][xx], (double) testdata, 1e-20, 1e-10));
            continue;
          }
          BOOST_REQUIRE(fcmp((double) image_data[zz][yy][xx], 0., 1e-20, 1e-10));
        }
      }
    }
    DMDAVecRestoreArray(warped->dmda(), warped->data_vector(), &image_data);
  }
}

BOOST_AUTO_TEST_CASE(shear_warp)
{

  for(integer curr_dof = 0; curr_dof < map->ndof(); curr_dof++)
  {
    std::cout << "Iteration " << curr_dof << "\n";
    // reset map
    map = MapBase::make_map_for_mask(*mask, nodespacing);
    image = std::make_unique<Image>(imgshape, MPI_COMM_WORLD);

    intcoord image_local_lo = image->local_offset();
    intcoord image_local_hi = image_local_lo + image->local_shape();
    // Access data via dmda and insert value
    floating*** image_data;
    DMDAVecGetArray(image->dmda(), image->data_vector(), &image_data);
    intcoord curr_loc;
    auto& xx = curr_loc[0];
    auto& yy = curr_loc[1];
    auto& zz = curr_loc[2];
    for(zz=image_local_lo[2]; zz < image_local_hi[2]; zz++)
    {
      for(yy=image_local_lo[1]; yy < image_local_hi[1]; yy++)
      {
        for(xx=image_local_lo[0]; xx < image_local_hi[0]; xx++)
        {
          if(curr_loc >= block_offset && curr_loc < (block_offset + block_shape))
          {
            image_data[zz][yy][xx] = testdata;
          }
        }
      }
    }
    DMDAVecRestoreArray(image->dmda(), image->data_vector(), &image_data);
    image->update_local_vector();

    // Set one dimension of the map to shear
    floating**** map_data;
    DMDAVecGetArrayDOF(map->dmda(), map->global_vector(), &map_data);
    intcoord map_local_lo = map->local_offset();
    intcoord map_local_hi = map_local_lo + map->local_shape();
    for(integer zz=map_local_lo[2]; zz<map_local_hi[2]; zz++)
    {
      for(integer yy=map_local_lo[1]; yy<map_local_hi[1]; yy++)
      {
        for(integer xx=map_local_lo[0]; xx<map_local_hi[0]; xx++)
        {
          map_data[zz][yy][xx][curr_dof] = xx * offset;
        }
      }
    }
    DMDAVecRestoreArrayDOF(map->dmda(), map->global_vector(), &map_data);

    auto warped = image->warp(*map);

    std::ostringstream fnamess;
    fnamess << "shear" << curr_dof << ".xdmf:/warped";
    auto wtr = XDMFWriter(fnamess.str(), MPI_COMM_WORLD);
    wtr.write_image(*warped);
  }
}


BOOST_AUTO_TEST_SUITE_END()
