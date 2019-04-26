#!/usr/bin/env python3

import argparse
import os
import sys

from .application_routines import run_pfire
from .analysis_routines import compare_image_results

def parse_args():
    """
    Grab minimal command line config
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("pfire_config")

    return parser.parse_args()


def main():
    if sys.version_info < (3, 5):
        print("This script requires minimum python version 3.5",
              file=sys.stderr)
        sys.exit(-1)

    args = parse_args()

    pfire_result = run_pfire(args.pfire_config)
    path_parts = os.path.split(pfire_result.registered_path)
    regression_image = os.path.join(path_parts[0],
                                    "accepted_result",
                                    path_parts[1])

    image_comparison = compare_image_results(pfire_result.fixed_path,
                                             pfire_result.moved_path,
                                             pfire_result.registered_path,
                                             regression_image)

if __name__ == "__main__":
    main()
