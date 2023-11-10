#!/usr/bin/env python3

import argparse
import os
import sys

from docutils.core import publish_string

from .analysis_routines import compare_image_results, compare_map_results
from .testinstance import TestInstance
from .application_routines import pFIRERunnerMixin, ShIRTRunnerMixin

def parse_args():
    """
    Grab minimal command line config
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("pfire_config")
    parser.add_argument("--test_name", type=str, default="analysis")

    return parser.parse_args()


def main():
    """ Run a single regression test
    """
    if sys.version_info < (3, 5):
        print("This script requires minimum python version 3.5",
              file=sys.stderr)
        sys.exit(-1)

    args = parse_args()

    test = ComparisonTest(args.pfire_config, name=args.test_name)

    test.run()
    res = test.generate_report()

    print("HTML report written to \"{}\"".format(res))


class ComparisonTest(TestInstance, pFIRERunnerMixin, ShIRTRunnerMixin):
    """ Test by comparing with results from ShIRT
    """

    def __init__(self, pfire_config, name=None, output_path=None):

        super().__init__(pfire_config, name=name, output_path=output_path)
        self.run_errstring = None
        print("output_path: {}".format(self.output_path))
        
        print("\n pFIRE config=" +pfire_config)
        exit

    def run(self):
        """ Run pfire against provided config
        """
        # First run pFIRE
        try:
            self.run_pfire(self.pfire_config)
        except RuntimeError as err:
            self.run_errstring = str(err)
            return False

        # Now run ShIRT
        try:
            self.run_shirt(self.pfire_config)
        except RuntimeError as err:
            self.run_errstring = str(err)
            return False

        return True


    def generate_report(self):
        """ Generate rst for report on ShiRT/pFIRE comparison
        """
        tables = ""
        images = ""

        tbl, img = compare_image_results(self.pfire_fixed_path,
                                         self.pfire_moved_path,
                                         self.shirt_reg_path,
                                         self.pfire_reg_path,
                                         fig_dir=self.output_path,
                                         cmpname="ShIRT")
        tables = "\n".join((tables, tbl))
        images = "\n".join((images, img))

        tbl, img = compare_map_results(self.shirt_map_path,
                                       self.pfire_map_path,
                                       fig_dir=self.output_path,
                                       cmpname="ShIRT")
        tables = "\n".join((tables, tbl))
        images = "\n".join((images, img))

        rst = "\n".join((tables, images))
        with open(self.report_file, 'wt') as fh:
            fh.write(publish_string(rst, writer_name='html5').decode())

        return self.report_file


if __name__ == "__main__":
    main()
