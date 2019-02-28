#define BOOST_TEST_MODULE shirtloader 
#include "test_common.hpp"

#include <algorithm>
#include <numeric>

#include "types.hpp"
#include "image.hpp"
#include "iterator_routines.hpp"

struct im
{
  im()
  {
    MPI_Comm_size(PETSC_COMM_WORLD, &commsize);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  }

  int rank;
  int commsize;
  std::string datapath = "shirtdata.image";
};

BOOST_FIXTURE_TEST_SUITE(shirtloader, im)

  BOOST_AUTO_TEST_CASE(test_natural_order_load)
  {
    // Load test data
    std::unique_ptr<Image> img = Image::load_file(datapath, nullptr, PETSC_COMM_WORLD);

    // img->global_vec() is petsc ordered, convert to natural:
    Vec_unique natvec(create_unique_vec());
    DMDACreateNaturalVector(*img->dmda(), natvec.get());
    DMDAGlobalToNaturalBegin(*img->dmda(), *img->global_vec(), INSERT_VALUES, *natvec);
    DMDAGlobalToNaturalEnd(*img->dmda(), *img->global_vec(), INSERT_VALUES, *natvec);

    // Natural ordered data should just run from 1 to N
    integer vec_lo, vec_hi, local_len;
    VecGetOwnershipRange(*natvec, &vec_lo, &vec_hi);
    local_len = vec_hi - vec_lo;

    floatvector cmp(local_len);
    std::iota(cmp.begin(), cmp.end(), vec_lo);
    
    floating *ptr;
    VecGetArray(*natvec, &ptr);

    BOOST_CHECK(std::equal(cmp.cbegin(), cmp.cend(), ptr)); 
  }

BOOST_AUTO_TEST_SUITE_END()
