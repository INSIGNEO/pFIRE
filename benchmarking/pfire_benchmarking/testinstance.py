#!/usr/bin/env python3

""" Abstract test classes
"""

import os

class TestInstance:
    """ Base class for specific test instances
    """

    def __init__(self, pfire_config, name=None, output_path=None):
        if name is not None:
            self.name = name
        else:
            self.name = os.path.splitext(os.path.basename(pfire_config))[0]

        if output_path is None:
            self.output_path = os.path.normpath('.')
        else:
            self.output_path = os.path.normpath(output_path)

        self.fig_dir = os.path.join(self.output_path, self.name)

        self.pfire_config = pfire_config

    def run(self):
        """ Template for test runner
        """
        raise NotImplementedError("testinstance should be subclassed")

    def generate_report(self):
        """ Template for test report generator
        """
        raise NotImplementedError("testinstance should be subclassed")
