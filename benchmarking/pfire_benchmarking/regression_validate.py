#!/usr/bin/env python3

import argparse
import sys

from .analysis_routines import compare_image_results
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

    return parser.parse_args()


def main():
    if sys.version_info < (3, 5):
        print("This script requires minimum python version 3.5",
              file=sys.stderr)
        sys.exit(-1)

    args = parse_args()

    test = ComparisonTest(args.pfire_config, accepted_image=args.accepted_image,
                          accepted_map=args.accepted_image)

    test.run()
    test.generate_report()


class ComparisonTest(TestInstance, pFIRERunnerMixin):
    """ Test by comparing with an accepted result image/map
    """

    def __init__(self, pfire_config, name=None, accepted_image=None,
                 accepted_map=None):
        # Check this before we go anywhere else
        if not (accepted_map or accepted_image):
            raise ValueError("At least one of accepted_image or accepted_map "
                             "must be provided")
        super().__init__(pfire_config, name=name)


        self.accepted_image_path = accepted_image
        self.accepted_map_path = accepted_map
        self.run_errstring = None

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
        
        image_comparison = compare_image_results(self.pfire_fixed_path,
                                                 self.pfire_moved_path,
                                                 self.pfire_reg_path,
                                                 self.accepted_image_path)




if __name__ == "__main__":
    main()
