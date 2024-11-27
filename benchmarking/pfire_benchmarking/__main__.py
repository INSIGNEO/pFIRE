#!/usr/bin/env python3

""" Testsuite entry point
"""

import argparse
import os

from .testdespatcher import TestDespatcher

def _parse_args():
    parser = argparse.ArgumentParser(description="Run pFIRE integration tests")

    parser.add_argument('dir', nargs='?', default=os.getcwd(), metavar="datadir",
                        help="Test data directory will be inspected recursively")

    parser.add_argument('--pfire-executable', nargs='?', default="pfire", metavar="pfire_exec_filename",
                        dest="pfire_exec_filename",
                        help="pFIRE executable filepath (filename with full path)")

    
    parser.add_argument('--output', '-o', default='.', metavar="outdir",
                        dest="output_dir",
                        help="Path at which to output results")

    return parser.parse_args()


def main():
    """ Run the testsuite over a directory tree
    """
    args = _parse_args()
    print(args)
    print("pfire_exec_filename=", args.pfire_exec_filename )
    
    testsuite = TestDespatcher(output_dir=args.output_dir, pfire_exec_filename=args.pfire_exec_filename ) 

    if not os.path.exists(args.dir):
        print("Error, path {} does not exist".format(args.dir))
        return
    # Add tests or tests in directory tree
    if os.path.isfile(args.dir):
        testsuite.add_test(args.dir)
    else:
        ntests = testsuite.find_tests(args.dir)
        if ntests == 0:
            print("No tests found.")
            return

    testsuite.run_tests()

    testsuite.create_aggregate_report()


if __name__ == "__main__":
    main()
