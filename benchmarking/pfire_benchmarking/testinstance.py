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
            self.output_path = os.path.normpath(os.path.join('.', self.name))
        else:
            self.output_path = os.path.normpath(os.path.join(output_path, self.name))

        os.makedirs(self.output_path, exist_ok=True)

        self.report_file = os.path.join(self.output_path, self.name + '.html')

        self.pfire_config = pfire_config

    def run(self):
        """ Template for test runner
        """
        raise NotImplementedError("testinstance should be subclassed")

    def generate_report(self):
        """ Template for test report generator
        """
        raise NotImplementedError("testinstance should be subclassed")
