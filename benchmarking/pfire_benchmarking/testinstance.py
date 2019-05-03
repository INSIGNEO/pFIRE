#!/usr/bin/env python3

import os

from configobj import ConfigObj

class TestInstance(object):
    """ Base class for specific test instances
    """

    def __init__(self, pfire_config, name=None):
        if name is not None:
            self.name = name
        else:
            self.name = os.path.splitext(os.path.basename(pfire_config))[0]

        self.pfire_config = pfire_config

    def run(self):
        raise NotImplementedError("testinstance should be subclassed")

    
    def generate_report(self):
        raise NotImplementedError("testinstance should be subclassed")




