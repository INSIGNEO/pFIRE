#define BOOST_TEST_MODULE map
#include "test_common.hpp"


#include "types.hpp"
#include "map.hpp"
#include "image.hpp"
#include "hdfwriter.hpp"

struct envobjs
{
  envobjs() : image(imgshape), map(image, nodespacing)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Create vec for addition values
    Vec &disps = map.displacements(); //values irrelevant but right shape so duplicate
    Vec_unique addvec(create_unique_vec());
    VecDuplicate(disps, addvec.get());

    // Now set values 0-n per dim
    if(rank == 0)
    {
      integer vecsize, mapsize;
      VecGetSize(*addvec, &vecsize);
      mapsize = map.size();
      intvector locs(vecsize);
      floatvector vals(vecsize);
      std::iota(locs.begin(), locs.end(), 0);
      std::transform(locs.begin(), locs.end(), vals.begin(),
                     [mapsize](integer loc) ->floating{return loc%mapsize;});
      VecSetValues(*addvec, vecsize, locs.data(), vals.data(), INSERT_VALUES);
    }
    VecAssemblyBegin(*addvec);
    VecAssemblyEnd(*addvec);

    // Map is initially zero, so an update should effectively set values
    map.update(*addvec); 
  }
  int rank;
  intvector imgshape = {8, 5, 1};
  floatvector nodespacing = {3, 3, 1};
  Image image;
  Map map;
};


BOOST_FIXTURE_TEST_SUITE(map, envobjs)

  BOOST_AUTO_TEST_CASE(test_map_get_single_dim_natural)
  {
    map.print_dm_info();
    Vec_unique natvec(map.get_single_dim_natural(0));
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_NATIVE);
    VecView(*natvec, PETSC_VIEWER_STDOUT_WORLD);
  }

  BOOST_AUTO_TEST_CASE(test_map_get_single_dim_petsc)
  {
    map.print_dm_info();
    Vec_unique petscvec(map.get_single_dim_petsc(0));
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_NATIVE);
    VecView(*petscvec, PETSC_VIEWER_STDOUT_WORLD);
  }

  BOOST_AUTO_TEST_CASE(test_map_get_raw_data_row_major)
  {
    Vec_unique natvec(map.get_raw_data_row_major(0));
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_NATIVE);
    VecView(*natvec, PETSC_VIEWER_STDOUT_WORLD);
    HDFWriter wtr("rmdata.h5", MPI_COMM_WORLD);
    wtr.write_map(map);
  }

BOOST_AUTO_TEST_SUITE_END()
