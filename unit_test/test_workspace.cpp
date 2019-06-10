#define BOOST_TEST_MODULE workspace
#include "test_common.hpp"

#include <petscdmda.h>

#include "image.hpp"
#include "map.hpp"
#include "types.hpp"
#include "workspace.hpp"

struct envobjs {
  envobjs()
    : image(imgshape), mask(Mask::full_image(image)), map(nodespacing, *mask), workspace(image, map)
  {
  }

  intvector imgshape = {5, 9, 7};
  floatvector nodespacing = {3, 3, 3};
  Image image;
  std::unique_ptr<Mask> mask;
  Map map;
  WorkSpace workspace;
};

BOOST_FIXTURE_TEST_SUITE(indexing, envobjs)

BOOST_AUTO_TEST_CASE(test_scatter_to_stacked)
{
  for (size_t idx = 0; idx < workspace.m_globaltmps.size(); idx++)
  {
    VecSet(*workspace.m_globaltmps[idx], (floating)idx);
    VecAssemblyBegin(*workspace.m_globaltmps[idx]);
    VecAssemblyEnd(*workspace.m_globaltmps[idx]);
  }
  workspace.scatter_grads_to_stacked();
  integer offset = 0;
  for (size_t idx = 0; idx < workspace.m_globaltmps.size(); idx++)
  {
    integer vecsize, lo, hi, llo, lhi;
    PetscErrorCode perr = VecGetSize(*workspace.m_globaltmps[idx], &vecsize);
    CHKERRXX(perr);
    perr = VecGetOwnershipRange(*workspace.m_globaltmps[idx], &llo, &lhi);
    CHKERRXX(perr);
    lo = offset;
    hi = offset + vecsize;
    if (lo >= lhi || hi < llo)
    {
      continue;
    }
    lo = lo > llo ? lo : llo;
    hi = hi < lhi ? hi : lhi;
    integer localsize = hi - lo;
    intvector idxn(localsize, lo);
    std::iota(idxn.begin(), idxn.end(), lo);
    floatvector data(localsize, 0);
    perr = VecGetValues(*workspace.m_stacktmp, localsize, idxn.data(), data.data());
    BOOST_REQUIRE(std::all_of(
        data.begin(), data.end(), [idx](floating x) -> bool { return x == floating(idx); }));
    offset += vecsize;
  }
}

BOOST_AUTO_TEST_SUITE_END()
