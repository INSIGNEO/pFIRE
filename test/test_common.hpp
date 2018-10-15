#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

#include <boost/test/unit_test.hpp>

#include<petscsys.h>

struct PetscSetup {
  PetscSetup() { PetscInitialize(nullptr, nullptr, nullptr, nullptr);}
  ~PetscSetup() { PetscFinalize();}
};

BOOST_GLOBAL_FIXTURE(PetscSetup);

#endif
