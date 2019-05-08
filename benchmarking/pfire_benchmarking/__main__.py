#!/usr/bin/env python3

""" Testsuite entry point
"""

import argparse
import os

from .testdespatcher import TestDespatcher

def _parse_args():
    parser = argparse.ArgumentParser(description="Run pFIRE integration tests")

    parser.add_argument('dir', nargs='?', default=os.getcwd(), metavar="datadir",
                        help="Test data directory, will be inspected recursively")
    parser.add_argument('--output', '-o', default='.', metavar="outdir",
                        help="Path at which to output results")

    return parser.parse_args()


def main():
    """ Run the testsuite over a directory tree
    """
    args = _parse_args()

    testsuite = TestDespatcher(output_dir=args.output)

    testsuite.find_tests(args.dir)

    testsuite.run_tests()

    testsuite.create_aggregate_report()


if __name__ == "__main__":
    main()
