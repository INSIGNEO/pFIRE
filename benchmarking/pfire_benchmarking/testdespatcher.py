#!/usr/bin/env python3

""" Despatcher class for collecting and running tests
"""

from pathlib import PurePosixPath
import os
from textwrap import dedent
from docutils.core import publish_string
from configobj import ConfigObj

from .regression_validate import RegressionTest
from .cross_validate import ComparisonTest

class TestDespatcher:
    """ Aggregates and runs a collection of tests
    """

    test_types = {"regression": RegressionTest,
                  "comparison": ComparisonTest}

    def __init__(self, output_dir=None):
        self.tests = []
        if output_dir:
            self.output_dir = os.path.normpath(output_dir)
        else:
            self.output_dir = os.path.normpath('.')


    def add_test(self, testconfig_path):
        """ Add a new test by parsing configfile
        """
        testdir = os.path.dirname(testconfig_path)
        # opening explicitly causes failure on file nonexistence
        with open(testconfig_path, 'r') as fh:
            testconfig = ConfigObj(fh)

        try:
            testiniter = self.test_types[testconfig['type'].lower()]
        except KeyError:
            types_str = ", ".join(self.test_types.keys())
            raise RuntimeError("Test type (\"type\")must be one of \"{}\""
                               "".format(types_str))

        testkwargs = {}
        if testconfig['type'].lower() == "regression":
            try:
                testkwargs['accepted_image'] = os.path.join(
                    testdir, testconfig['accepted_image'])
            except KeyError:
                pass
            try:
                testkwargs['accepted_map'] = os.path.join(
                    testdir, testconfig['accepted_map'])
            except KeyError:
                pass
            if not testkwargs:
                raise RuntimeError("Regression test must specify at least one "
                                   "accepted result file")

        try:
            testkwargs['name'] = testconfig['name']
        except KeyError:
            testkwargs['name'] = None

        testkwargs['output_path'] = self.output_dir

        pfire_config_path = os.path.join(testdir, testconfig['pfire_config'])
        try:
            test = testiniter(pfire_config_path, **testkwargs)
        except KeyError:
            raise RuntimeError("A pFIRE configuration file (\"pfire_config\") "
                               "must be specified")

        print("Found \"{}\"".format(test.name))
        self.tests.append(test)


    def find_tests(self, search_dir):
        """ Find all tests in a directory tree
        """
        for dirname, _, fnames in os.walk(search_dir):
            for fname in fnames:
                if fname.endswith(".testconf"):
                    testpath = os.path.join(dirname, fname)
                    try:
                        self.add_test(testpath)
                    except RuntimeError as err:
                        print("Error adding {}: {}".format(testpath, err))

    def run_tests(self):
        """ Run all tests producing report files
        """
        for test in self.tests:
            test.run()
            test.generate_report()


    def create_aggregate_report(self):
        """ Create summary report page linking to individual tests
        """
        index_rst = []

        index_rst.append(dedent("""\
            =================
            Testsuite Results
            =================
            """))

        for test in self.tests:
            report_relpath = PurePosixPath(test.report_file).relative_to(self.output_dir)
            index_rst.append(dedent("""\
                `{0}`_

                .. _{0}: {1}
            """.format(test.name, report_relpath)))

        index_rst = "\n".join(index_rst)

        index_html = publish_string(index_rst, writer_name='html5')

        with open(os.path.join(self.output_dir, "report.html"), 'wt') as fh:
            fh.write(index_html.decode())
