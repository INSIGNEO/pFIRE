#!/usr/bin/env python3

import argparse
import sys

from .application_routines import run_pfire, run_shirt
from .analysis_routines import compare_map_results, compare_image_results

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
    shirt_result = run_shirt(args.pfire_config)

    image_comparison = compare_image_results(shirt_result.fixed_path,
                                             shirt_result.moved_path,
                                             shirt_result.registered_path,
                                             pfire_result.registered_path)

    map_coeffs = compare_map_results(shirt_result.map_path,
                                     pfire_result.map_path,
                                     save_figs=args.pfire_config)

    print(image_comparison)
    print(map_coeffs)


if __name__ == "__main__":
    main()
