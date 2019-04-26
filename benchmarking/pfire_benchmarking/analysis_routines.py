#!/usr/bin/env python3

from textwrap import wrap
import os

import matplotlib.pyplot as plt
import scipy.stats as sps

import flannel.io as fio

from .image_routines import (calculate_mutual_information,
                             load_image,
                             load_pfire_map)

def compare_image_results(fixed_path, moved_path, shirt_registered_path,
                          pfire_registered_path, save_figs=None):
    """Compare ShIRT and pFIRE registered images
    """

    fixed_data = load_image(fixed_path)
    moved_data = load_image(moved_path)
    shirt_data = load_image(shirt_registered_path)
    pfire_data = load_image(pfire_registered_path)

    mi_start = calculate_mutual_information(fixed_data, moved_data,
                                            return_hist=True)
    mi_shirt = calculate_mutual_information(fixed_data, shirt_data,
                                            return_hist=True)
    mi_pfire = calculate_mutual_information(fixed_data, pfire_data,
                                            return_hist=True)
    mi_comparison = calculate_mutual_information(shirt_data, pfire_data,
                                                 return_hist=True)

    prof_start = mi_start[0]/mi_start[1]
    prof_shirt = mi_shirt[0]/min(mi_shirt[1], mi_shirt[2])
    prof_pfire = mi_pfire[0]/min(mi_pfire[1], mi_pfire[2])
    rdcy_norm = mi_comparison[0]/min(mi_comparison[1], mi_comparison[2])


    print("Normalized mutual information (proficiency):")
    print("Pre registration: {:.3f}".format(prof_start))
    print("Shirt registration: {:.3f}".format(prof_shirt))
    print("pFIRE registration: {:.3f}".format(prof_pfire))
    print("Similarity between results (normalized redundancy):")
    print("ShIRT vs. pFIRE: {:.3f}\n".format(rdcy_norm))

    if save_figs:
        plt.matshow(mi_start[-1], origin='lower')
        plt.title("\n".join(wrap("Pre-registration normalized mutual information: "
                                 "{:0.3f}".format(prof_start), 40)))
        plt.savefig("{}_prereg.png".format(os.path.splitext(save_figs)[0]))
        plt.close()

        plt.matshow(mi_shirt[-1], origin='lower')
        plt.title("\n".join(wrap("ShIRT registration normalized mutual information: "
                                 "{:0.3f}".format(prof_shirt), 40)))
        plt.savefig("{}_shirt.png".format(os.path.splitext(save_figs)[0]))
        plt.close()

        plt.matshow(mi_pfire[-1], origin='lower')
        plt.title("\n".join(wrap("pFIRE registration normalized mutual information: "
                                 "{:0.3f}".format(prof_pfire), 40)))
        plt.savefig("{}_pfire.png".format(os.path.splitext(save_figs)[0]))
        plt.close()

        plt.matshow(mi_comparison[-1], origin='lower')
        plt.title("\n".join(wrap("ShIRT - pFIRE comparison (normalized redundancy): "
                                 "{:0.3f}".format(rdcy_norm), 40)))
        plt.savefig("{}_comparison.png".format(os.path.splitext(save_figs)[0]))
        plt.close()

    mi_data = (prof_start, prof_shirt, prof_pfire, rdcy_norm)

    return mi_data


def compare_map_results(shirt_map_path, pfire_map_path, save_figs=None):
    """Compare ShIRT and pFIRE displacement maps
    """

    shirt_map = fio.load_map(shirt_map_path)[0][0:3]
    pfire_map = load_pfire_map(pfire_map_path)

    print("Map coefficients of determination (R^2), per dimension:")

    corr_data = []

    plt.plot(shirt_map[0].flatten(), marker='x', ls='none', label="ShIRT")
    plt.plot(pfire_map[0].flatten(), marker='+', ls='none', label="pFIRE")
    corr_data.append(sps.linregress(shirt_map[0].flatten(),
                                    pfire_map[0].flatten())[2])
    print("X: {:0.3}".format(corr_data[-1]**2))
    plt.title("Map X component, R^2={:0.3}".format(corr_data[-1]**2))
    plt.legend()
    if save_figs:
        plt.savefig("{}_map_x.png".format(os.path.splitext(save_figs)[0]))
    plt.close()

    plt.plot(shirt_map[1].flatten(), marker='x', ls='none', label="ShIRT")
    plt.plot(pfire_map[1].flatten(), marker='+', ls='none', label="pFIRE")
    corr_data.append(sps.linregress(shirt_map[1].flatten(),
                                    pfire_map[1].flatten())[2])
    plt.title("Map Y component, R^2={:0.3}".format(corr_data[-1]**2))
    print("Y: {:0.3}".format(corr_data[-1]**2))
    plt.legend()
    if save_figs:
        plt.savefig("{}_map_y.png".format(os.path.splitext(save_figs)[0]))
    plt.close()

    try:
        plt.plot(shirt_map[2].flatten(), marker='x', ls='none', label="ShIRT")
        plt.plot(pfire_map[2].flatten(), marker='+', ls='none', label="pFIRE")
        corr_data.append(sps.linregress(shirt_map[2].flatten(),
                                        pfire_map[2].flatten())[0])
        plt.title("Map Z component, R^2={:0.3}".format(corr_data[2]**2))
        print("Z: {:0.3}".format(corr_data[2]**2))
        plt.legend()
        if save_figs:
            plt.savefig("{}_map_z.png".format(os.path.splitext(save_figs)[0]))
        plt.close()
    except IndexError:
        pass

    return corr_data
