#!/usr/bin/env python3


import os
from textwrap import dedent
from docutils.core import publish_string
from configobj import ConfigObj

from . import RegressionTest, ComparisonTest

class TestDespatcher:
    """ Aggregates and runs a collection of tests
    """

    test_types = {"regression": RegressionTest,
                  "comparison": ComparisonTest}

    def __init__(self):
        self.tests = []

    def add_test(self, testconfig):
        """ Add a new test by parsing configfile
        """
        with open(testconfig, 'r') as fh:
            testconfig = ConfigObj(fh)

        try:
            test = self.test_types[testconfig['type'].lower()](testconfig)
        except KeyError:
            types_str = ", ".join(self.test_types.keys())
            raise ValueError("Test type must be one of \"{}\""
                             "".format(types_str))

        self.tests.append(test)


    def find_tests(self, search_dir):
        """ Find all tests in a directory tree
        """
        for dirname, _, fnames in os.walk(search_dir):
            for fname in fnames:
                if fname.endswith(".testconf"):
                    self.add_test(os.path.join(dirname, fname))

    def run_tests(self):
        """ Run all tests producing report files
        """
        for test in self.tests:
            test.run()
            test.generate_report()


    def create_aggregate_report(self):
        """ Create summary report page linking to individual tests
        """
        result_files = {}
        for test in self.tests:
            result_files[test.name] = test.report_filename

        index_rst = []

        index_rst.append(dedent("""\
            =================
            Testsuite Results
            =================
            """))

        for test_name, result_file in result_files.items():
            index_rst.append(dedent("""\
                `{0}`_

                .. _{0}: {1}
            """.format(test_name, result_file)))

        index_rst = "\n".join(index_rst)

        with open("index.rst", 'wt') as fh:
            fh.write(index_rst)

        index_html = publish_string(index_rst, writer_name='html5')

        with open("index.html", 'wt') as fh:
            fh.write(index_html.decode())
