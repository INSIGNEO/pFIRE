#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

#include <boost/test/unit_test.hpp>

#include "setup.hpp"

struct pFIRESetup {
  pFIRESetup() { pfire_setup({}, true);}
  ~pFIRESetup() { pfire_teardown();}
};

BOOST_GLOBAL_FIXTURE(pFIRESetup);

#endif
