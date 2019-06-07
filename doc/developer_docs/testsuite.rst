================================
Using the Integration Test Suite
================================

Along with the pFIRE library and executable, the codebase also includes unit tests and integration
tests for validation, benchmarking and validation of pFIRE's functionality.  The use of these is
documented in detail here.


Integration Testing
===================

Integration testing involves running the complete program against various input data and comparing
either to an accepted result, either in the form of predetermined result file or the output of
another program (this could be an earlier version of pFIRE or a different image registration code).

Integration testing for pFIRE can be performed using the `pfire_benchmarking` python module,
located in the `benchmarking/` subdirectory of the source code repository.


Installation
------------

Since it is provided as
a standard python module the test suite may be installed on any system with python version 3.5 or
later, either using pip or by manually running the provided setup.py:

.. code-block:: console

  $ cd benchmarking
  # then
  $ pip install --user .
  #or
  $ python3 setup.py

The test suite is dependent on several common python packages. If using pip these should be
installed for you automatically.  If installing manually it is up to you to install them.

Once installed, the testsuite may be run using the `pfire-integration-test` command.


Running Integration Tests
-------------------------

The integration test suite can be run either on a single test file or on all tests in a directory
tree. There are two types of test that the suite supports, comparison with a result file, or by
running both pFIRE and its predecessor (ShIRT) against the same test data and comparing the
results.

To run against a single configuration file, pass that test file as a parameter to the testsuite

.. code-block:: console

  $ pfire-integration-test /path/to/testfile.testconf

or to run against all test files in a directory tree, pass the root directory of that tree

.. code-block:: console

  $ pfire-integration-test /path/to/root/of/testtree/

the testsuite will then be run against the supplied test(s), and an html report generated detailing
the results.  Note that for a test configuration file to be detected in the test tree it's name
should end with `.testconf`.

Included Tests
--------------

A set of comparison and regression tests are included with the pFIRE source code, in the
`testdata/integration` subdirectory.

Integration Test Configuration
------------------------------

Integration tests are described and configured using a short configuration file, with `ini` syntax.
This configuration file should have the extension `.testconf` in order to be detected by the
testsuite when part of a tree of tests.

The complete syntax of the config file is:

.. code-block:: ini

  name = # string: Name of the Test
  type = # string: one of "comparison" or "regression"
  pfire_config = # string: the pfire config file to use
  accepted_image = # string: the image to compare with ("regression" mode only)
  accepted_map = # string: the map to compare with ("regression" mode only)

The test suite will run pFIRE using the supplied configuration file, and compare the output to
either the provided file(s) (`type=regression`) or to the results of running ShIRT with the same
configuration (`type=comparison`).  If the type is specified as `comparison` then the `ShIRT`
binary must be located on the system `$PATH`.

