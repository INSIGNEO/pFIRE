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

#define BOOST_TEST_MODULE hdfwriter
#include "test_common.hpp"
#include "test_helpers.hpp"

#include <petscdmda.h>

#include "image.hpp"
#include "mapbase.hpp"
#include "mask.hpp"
#include "types.hpp"
#include "hdfwriter.hpp"

class envobjs {
public:
  envobjs()
    : image(std::make_unique<Image>(imgshape, MPI_COMM_WORLD)), mask(Mask::create_filled_mask(*image)),
      map(MapBase::make_map_for_mask(*mask, nodespacing))
  {
  }

  const floating testdata = 1.0;
  const integer offset = 5;
  const intcoord imgshape = {150, 100, 50};
  const intcoord nodespacing = {30, 30, 30};
  const intcoord block_shape = {30,30,30};
  const intcoord block_offset = {10,10,10};
  std::shared_ptr<Image> image;
  std::shared_ptr<Mask> mask;
  std::shared_ptr<MapBase> map;
};

BOOST_FIXTURE_TEST_SUITE(indexing, envobjs)

BOOST_AUTO_TEST_CASE(test_write_hdf)
{

  // Write dof number to all nodes in dof
  floating**** map_data;
  DMDAVecGetArrayDOF(map->dmda(), map->global_vector(), &map_data);
  intcoord map_local_lo = map->local_offset();
  intcoord map_local_hi = map_local_lo + map->local_shape();
  for(integer idof=0; idof<4; idof++)
    {
      for(integer zz=map_local_lo[2]; zz<map_local_hi[2]; zz++)
      {
        for(integer yy=map_local_lo[1]; yy<map_local_hi[1]; yy++)
        {
          for(integer xx=map_local_lo[0]; xx<map_local_hi[0]; xx++)
          {
            map_data[zz][yy][xx][idof] = idof;
          }
        }
      }
    }
  DMDAVecRestoreArrayDOF(map->dmda(), map->global_vector(), &map_data);

  {
    auto wtr = HDFWriter("testdata.h5:/map", MPI_COMM_WORLD); 
    wtr.write_map(*map);
  }

  // Now open and check HDF file, each dof should be all that dof number
  // No need to MPI this because we are doing read-only access
  hid_t file_h = H5Fopen("testdata.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  for(integer idof=0; idof < 3; idof++)
  {
    coord<hsize_t> local_shape = map->local_shape();
    coord<hsize_t> local_offset = map->local_offset();

    std::string dset_name = "/map/" + BaseWriter::components()[idof];
    hid_t dset_h = H5Dopen(file_h, dset_name.c_str(), H5P_DEFAULT);
    hid_t fspace_h = H5Dget_space(dset_h);
    H5Sselect_hyperslab(fspace_h, H5S_SELECT_SET, local_offset.data(), nullptr, local_shape.data(),
                        nullptr);
  
    hid_t mspace_h = H5Screate_simple(3, local_shape.data(), nullptr);

    auto* dataspace = new floating[map->num_owned_nodes()];
    H5Dread(dset_h, HDFWriter::get_hdf5_petsc_scalar(), mspace_h, fspace_h, H5P_DEFAULT, dataspace);

    for(integer iidx=0; iidx < map->num_owned_nodes(); iidx++)
    {
      BOOST_REQUIRE(fcmp((double) dataspace[iidx], double(idof+1), 1e-20, 1e-20));
    }

    delete[] dataspace;
  }

  H5Fclose(file_h);
}

BOOST_AUTO_TEST_SUITE_END()
