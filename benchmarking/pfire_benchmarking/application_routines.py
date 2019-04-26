#!/usr/bin/env python3

import os
import subprocess as sp

import numpy as np

from configobj import ConfigObj

import flannel.io as fio

default_mask_name = "default_mask.mask"
config_defaults = {"mask": None}

class ResultObject:
    """
    Small object to hold registration result info
    """

    def __init__(self, registered_path, map_path, log_path, fixed_path=None,
                 moved_path=None):
        self.registered_path = registered_path
        self.map_path = map_path
        self.log_path = log_path
        self.fixed_path = fixed_path
        self.moved_path = moved_path


def build_shirt_config(config_file):
    """
    Build ShIRT config object given pfire config file
    """
    # read pfire config and merge with default options to get full option list
    defaults = ConfigObj(config_defaults)
    config = ConfigObj(config_file)
    defaults.merge(config)
    return defaults

def run_pfire(config_file, comm_size=1):
    """
    Run pFIRE using provided config file
    """
    config = ConfigObj(config_file)
    print("Running pFIRE on {}".format(config_file))
    pfire_args = ['pfire', config_file]

    logfile_name = "{}_pfire.log".format(os.path.splitext(config_file)[0])
    with open(logfile_name, 'w') as logfile:
        res = sp.run(pfire_args, stdout=logfile, stderr=logfile)

    if res.returncode != 0:
        raise RuntimeError("Failed to run pFIRE, check log for details: {}"
                           "".format(logfile_name))

    reg_path = str()
    map_path = str()
    with open(logfile_name, 'r') as logfile:
        for line in logfile:
            if line.startswith("Saving registered image to "):
                reg_path = line.replace("Saving registered image to ", "").strip()
            if line.startswith("Saving map to "):
                map_path = line.replace("Saving map to ", "").strip()
    if not reg_path:
        raise RuntimeError("Failed to extract registered image path from log")
    if not map_path:
        raise RuntimeError("Failed to extract map path from log")

    return ResultObject(reg_path, map_path, logfile_name, config['fixed'],
                        config['moved'])


def run_shirt(config_file):
    """
    Run ShIRT using a pfire config file for input
    """
    print("Running ShIRT on {}".format(config_file))
    config = build_shirt_config(config_file)

    map_path = 'map.map'
    reg_path = 'registered.image'

    if config['mask'] is None:
        maskname = default_mask_name
        data = fio.load_image(config['fixed'])
        data = np.full(data.shape, 1.0)
        fio.save_image(data, maskname)
        config['mask'] = maskname

    shirt_args = ['ShIRT', 'Register', 'verbose',
                  'Fixed', config['fixed'].replace('.image', ''),
                  'Moved', config['moved'].replace('.image', ''),
                  'Mask', config['mask'].replace('.mask', ''),
                  'NodeSpacing', config['nodespacing'],
                  'Registered', reg_path.replace('.image', ''),
                  'Map', map_path.replace('.map', '')]

    logfile_name = "{}_shirt.log".format(os.path.splitext(config_file)[0])
    with open(logfile_name, 'w') as logfile:
        sp.run(['ShIRT', 'setpath', 'DataPath', '.'])
        res = sp.run(shirt_args, stdout=logfile, stderr=logfile)

    if res.returncode != 0:
        raise RuntimeError("Failed to run ShIRT, check log for details: {}"
                           "".format(logfile_name))

    return ResultObject(reg_path, map_path, logfile_name, config['fixed'],
                        config['moved'])
