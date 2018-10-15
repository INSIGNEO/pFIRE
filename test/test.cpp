#define BOOST_TEST_MODULE MyCMakeTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test1) {
  BOOST_CHECK_EQUAL(42, 42);
}

BOOST_AUTO_TEST_CASE(test2) {
  BOOST_CHECK_EQUAL(42, 42);
}

BOOST_AUTO_TEST_CASE(test3) {
  BOOST_CHECK_EQUAL(42, 41);
}
