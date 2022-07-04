#!/usr/bin/env python3

import argparse
import os
import sys

from docutils.core import publish_string

from .analysis_routines import compare_image_results, compare_map_results
from .testinstance import TestInstance
from .application_routines import pFIRERunnerMixin

def parse_args():
    """
    Grab minimal command line config
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("pfire_config")
    parser.add_argument("--accepted_image", type=str)
    parser.add_argument("--accepted_map", type=str)
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

    test = RegressionTest(args.pfire_config, 
                          name=args.test_name,
                          accepted_image=args.accepted_image,
                          accepted_map=args.accepted_map)

    test.run()
    res = test.generate_report()

    print("HTML report written to \"{}\"".format(res))


class RegressionTest(TestInstance, pFIRERunnerMixin):
    """ Test by comparing with an accepted result image/map
    """

    def __init__(self, pfire_config, name=None, accepted_image=None,
                 accepted_map=None, output_path=None):
        
        # Check this before we go anywhere else
        if not (accepted_map or accepted_image):
            raise ValueError("At least one of accepted_image or accepted_map "
                             "must be provided")
        super().__init__(pfire_config, name=name, output_path=output_path)


        self.accepted_image_path = accepted_image
        self.accepted_map_path   = accepted_map
        self.run_errstring = None

        print("\n pFIRE config=" +pfire_config)


    def run(self):
        """ Run pfire against provided config
        """
        try:
            self.run_pfire(self.pfire_config)
        except RuntimeError as err:
            self.run_errstring = str(err)
            return False

        return True


    def generate_report(self):

        tables = ""
        images = ""

        if self.accepted_image_path:
            tbl, img = compare_image_results(self.pfire_fixed_path,
                                             self.pfire_moved_path,
                                             self.accepted_image_path,
                                             self.pfire_reg_path,
                                             fig_dir=self.output_path,
                                             cmpname="Accepted")
            tables = "\n".join((tables, tbl))
            images = "\n".join((images, img))

        if self.accepted_map_path:
            tbl, img = compare_map_results(self.accepted_map_path,
                                           self.pfire_map_path,
                                           fig_dir=self.output_path,
                                           cmpname="Accepted")
            tables = "\n".join((tables, tbl))
            images = "\n".join((images, img))

        rst = "\n".join((tables, images))
        with open(self.report_file, 'wt') as fh:
            fh.write(publish_string(rst, writer_name='html5').decode())

        return self.report_file


if __name__ == "__main__":
    main()
