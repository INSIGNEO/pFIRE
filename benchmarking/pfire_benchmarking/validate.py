#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys
import tempfile
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sps
import flannel.io as fio

from configobj import ConfigObj

from .image_routines import (calculate_mutual_information,
                             load_pfire_image,
                             load_pfire_map)

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

def parse_args():
    """
    Grab minimal command line config
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("pfire_config")

    return parser.parse_args()



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
    print ("Running pFIRE on {}".format(config_file))
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
                reg_path = line.replace("Saving registered image to ", "")
            if line.startswith("Saving map to "):
                map_path = line.replace("Saving map to ", "")
    if not reg_path:
        raise RuntimeError("Failed to extract registered image path from log")
    if not map_path:
        raise RuntimeError("Failed to extract map path from log")

    return ResultObject(reg_path, map_path, logfile_name)

def run_shirt(config_file):
    """
    Run ShIRT using a pfire config file for input
    """
    print ("Running ShIRT on {}".format(config_file))
    config = build_shirt_config(config_file)

    map_path = 'map.map'
    reg_path = 'registered.image'

    if config['mask'] is None:
        maskname = default_mask_name
        data = fio.load_image(config['fixed'])
        data = np.full(data.shape, 1.0)
        fio.save_image(data, maskname)
        config['mask'] = maskname

    shirt_args = ['ShIRT', 'Register',
                  'Fixed', config['fixed'].replace('.image', ''),
                  'Moved', config['moved'].replace('.image', ''),
                  'Mask', config['mask'].replace('.mask',''),
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

def main():

    if sys.version_info < (3, 5):
        print("This script requires minimum python version 3.5",
              file=sys.stderr)
        sys.exit(-1)

    args = parse_args()

    pfire_result = run_pfire(args.pfire_config)
    shirt_result = run_shirt(args.pfire_config)

    fixed_data = fio.load_image(shirt_result.fixed_path)/255
    moved_data = fio.load_image(shirt_result.moved_path)/255
    shirt_data = fio.load_image(shirt_result.registered_path)/255
    pfire_data = load_pfire_image(pfire_result.registered_path)

    mi_start = calculate_mutual_information(fixed_data, moved_data,
                                            return_hist=True)
    mi_shirt = calculate_mutual_information(fixed_data, shirt_data,
                                            return_hist=True)
    mi_pfire = calculate_mutual_information(fixed_data, pfire_data,
                                            return_hist=True)
    mi_comparison = calculate_mutual_information(shirt_data, pfire_data,
                                                 return_hist=True)

    prof_start = mi_start[0]/mi_start[1]
    prof_shirt = mi_shirt[0]/mi_shirt[1]
    prof_pfire = mi_pfire[0]/mi_pfire[1]
    rdcy_comparison = mi_comparison[0]/(mi_comparison[1]+mi_comparison[2])
    rdcy_max = min(mi_comparison[1], mi_comparison[2])/(mi_comparison[1]+mi_comparison[2])
    rdcy_norm = rdcy_comparison/rdcy_max


    print("Normalized mutual information (proficiency):")
    print("Pre registration: {:.3f}".format(prof_start))
    print("Shirt registration: {:.3f}".format(prof_shirt))
    print("pFIRE registration: {:.3f}".format(prof_pfire))
    print("Similarity between results (normalized redundancy):")
    print("ShIRT vs. pFIRE: {:.3f}\n".format(rdcy_norm))

    plt.matshow(mi_start[-1], origin='lower')
    plt.title("\n".join(wrap("Pre-registration normalized mutual information: "
                             "{:0.3f}".format(prof_start), 40)))
    plt.savefig("{}_prereg.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()

    plt.matshow(mi_shirt[-1], origin='lower')
    plt.title("\n".join(wrap("ShIRT registration normalized mutual information: "
                             "{:0.3f}".format(prof_shirt), 40)))
    plt.savefig("{}_shirt.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()

    plt.matshow(mi_pfire[-1], origin='lower')
    plt.title("\n".join(wrap("pFIRE registration normalized mutual information: "
                             "{:0.3f}".format(prof_pfire), 40)))
    plt.savefig("{}_pfire.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()

    plt.matshow(mi_comparison[-1], origin='lower')
    plt.title("\n".join(wrap("ShIRT - pFIRE comparison (normalized redundancy): "
                             "{:0.3f}".format(rdcy_norm), 40)))
    plt.savefig("{}_comparison.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()



    shirt_map = fio.load_map(shirt_result.map_path)[0][0:3]
    pfire_map = load_pfire_map(pfire_result.map_path)

    print("Map correlation coefficients, per dimension:")


    plt.plot(shirt_map[0].flatten(), marker='x', ls='none', label="ShIRT")
    plt.plot(pfire_map[0].flatten(), marker='+', ls='none', label="pFIRE")
    corr_info = sps.linregress(shirt_map[0].flatten(), pfire_map[0].flatten())
    print("X: {:0.3}".format(corr_info[2]))
    plt.title("Map X component, R={:0.3}".format(corr_info[2]))
    plt.legend()
    plt.savefig("{}_map_x.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()

    plt.plot(shirt_map[1].flatten(), marker='x', ls='none', label="ShIRT")
    plt.plot(pfire_map[1].flatten(), marker='+', ls='none', label="pFIRE")
    corr_info = sps.linregress(shirt_map[1].flatten(), pfire_map[1].flatten())
    plt.title("Map Y component, R={:0.3}".format(corr_info[2]))
    print("Y: {:0.3}".format(corr_info[2]))
    plt.legend()
    plt.savefig("{}_map_y.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()

    plt.plot(shirt_map[2].flatten(), marker='x', ls='none', label="ShIRT")
    plt.plot(pfire_map[2].flatten(), marker='+', ls='none', label="pFIRE")
    corr_info = sps.linregress(shirt_map[2].flatten(), pfire_map[2].flatten())
    plt.title("Map Z component, R={:0.3}".format(corr_info[2]))
    print("Z: {:0.3}".format(corr_info[2]))
    plt.legend()
    plt.savefig("{}_map_z.png".format(os.path.splitext(args.pfire_config)[0]))
    plt.close()

if __name__ == "__main__":
    main()
