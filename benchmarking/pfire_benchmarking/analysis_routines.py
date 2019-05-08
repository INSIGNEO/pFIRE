#!/usr/bin/env python3

""" Mathematical analysis functions for image and map comparison
"""

from collections import namedtuple
from textwrap import wrap
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps

from tabulate import tabulate

from .image_routines import load_image, load_map

MIResult = namedtuple("mi_result", ['mi', 'hist'])

def calculate_entropy(prob_dist):
    r"""
    Calculate Shannon entropy of the provided probability distribution

    Shannon Entropy is defined as
    $H(X) = \sum_n p_X(x)\log_2{p_X(x)}$
    """
    # First disregard all values where p_X == 0 to avoid nans from log(p_X)
    normed_prob_dist = prob_dist[prob_dist > 0]
    normed_prob_dist /= normed_prob_dist.sum()
    entropy = -np.sum(normed_prob_dist * np.log2(normed_prob_dist))

    return entropy


def calculate_mutual_information(data1, data2, resolution=50,
                                 return_hist=False):
    r"""
    Calculate mutual information using Shannon entropy of provided data.

    Mutual Information is defined as:
    MI(X, Y) = H(X) + H(Y) - H(X,Y)

    Where H(X), H(Y) and H(X,Y) are the Shannon entropies of the probabilities
    and the joint probability of the data.

    N.B it is assumed that the two datasets are independent.

    Returns a tuple of MI(X,Y), H(X), H(Y), H(X,Y)
    """
    jointmax = max(data1.max(), data2.max())
    # First calculate probability density
    bin_edges = np.linspace(0, 1, num=resolution)
    prob1_2, _, _ = np.histogram2d(data1.flatten()/jointmax,
                                   data2.flatten()/jointmax,
                                   bins=bin_edges, density=True)
    prob1 = np.sum(prob1_2, axis=1)
    prob2 = np.sum(prob1_2, axis=0)

    entropy1 = calculate_entropy(prob1)
    entropy2 = calculate_entropy(prob2)
    entropy1_2 = calculate_entropy(prob1_2)

    mutual_information = entropy1 + entropy2 - entropy1_2

    if return_hist:
        return (mutual_information, entropy1, entropy2, entropy1_2, prob1_2)
    else:
        return (mutual_information, entropy1, entropy2, entropy1_2)


def plot_2dhist(data, path, title):
    """ Helper function to plot 2d histogram and return rst inclusion command.
    """
    plt.matshow(data, origin='lower', cmap='gray')
    plt.title("\n".join(wrap(title, 40)))
    plt.savefig(path)
    plt.close()

    return ".. image:: {}\n".format(os.path.basename(path))


def calculate_proficiency(alpha, beta):
    """ Calculate proficiency (normalized mutual information) of an image pair
    """
    alpha_data = load_image(alpha)
    beta_data = load_image(beta)

    res = calculate_mutual_information(alpha_data, beta_data, return_hist=True)

    prof = res[0]/min(res[1], res[2])

    return MIResult(prof, res[-1])


def compare_image_results(fixed_path, moved_path, accepted_path,
                          pfire_path, fig_dir=None, cmpname="accepted"):
    """Compare ShIRT and pFIRE registered images
    """
    if fig_dir:
        os.makedirs(os.path.normpath(fig_dir), mode=0o755, exist_ok=True)
    else:
        fig_dir = os.path.normpath('.')

    mi_start = calculate_proficiency(fixed_path, moved_path)
    mi_accepted = calculate_proficiency(fixed_path, accepted_path)
    mi_pfire = calculate_proficiency(fixed_path, pfire_path)
    mi_comparison = calculate_proficiency(accepted_path, pfire_path)

    res_table = [("Normalized mutual information (proficiency):", ""),
                 ("Fixed vs. Moved:", "{:.3f}".format(mi_start.mi)),
                 ("{} vs. Fixed:".format(cmpname), "{:.3f}".format(mi_accepted.mi)),
                 ("pFIRE vs. Fixed:", "{:.3f}".format(mi_pfire.mi)),
                 ("pFIRE vs. {}:".format(cmpname), "{:.3f}\n".format(mi_comparison.mi))]

    print(tabulate(res_table, headers="firstrow", tablefmt='grid') + "\n")

    rst_output = []
    rst_output.append(tabulate(res_table, headers="firstrow", tablefmt="rst"))
    rst_output.append("") # table must be followed by blank line

    image_rst = []

    if fig_dir:
        image_rst.append(plot_2dhist(
            mi_start.hist, os.path.join(fig_dir, "prereg.png"),
            "Fixed vs. Moved normalized mutual information: "
            "{:0.3f}".format(mi_start.mi)))

        image_rst.append(plot_2dhist(
            mi_accepted.hist, os.path.join(fig_dir, "accepted.png"),
            "{} vs. Fixed normalized mutual information: "
            "{:0.3f}".format(cmpname, mi_accepted.mi)))

        image_rst.append(plot_2dhist(
            mi_pfire.hist, os.path.join(fig_dir, "pfire.png"),
            "pFIRE vs Fixed normalized mutual information: "
            "{:0.3f}".format(mi_pfire.mi)))

        image_rst.append(plot_2dhist(
            mi_comparison.hist, os.path.join(fig_dir, "comparison.png"),
            "pFIRE vs. {} normalized mutual information: "
            "{:0.3f}".format(cmpname, mi_comparison.mi)))

    return ("\n".join(rst_output), "\n".join(image_rst))


def compare_map_results(cmp_map_path, pfire_map_path, fig_dir=None,
                        cmpname='Accepted'):
    """Compare ShIRT and pFIRE displacement maps
    """
    if fig_dir:
        os.makedirs(os.path.normpath(fig_dir), mode=0o755, exist_ok=True)

    cmp_map = load_map(cmp_map_path)
    pfire_map = load_map(pfire_map_path)

    table_entries = [("Map coefficients of determination (R^2), by dimension:", "")]
    image_entries = []

    for didx, dim in enumerate(['X', 'Y', 'Z']):
        try:
            corr = sps.linregress(cmp_map[didx].flatten(),
                                  pfire_map[didx].flatten())[2]
            table_entries.append(("{}:".format(dim), "{:0.3}".format(corr**2)))

            if fig_dir:
                savepath = os.path.join(fig_dir, "map_{}.png".format(dim.lower()))
                plt.plot(cmp_map[didx].flatten(), marker='x', ls='none',
                         label=cmpname)
                plt.plot(pfire_map[didx].flatten(), marker='+', ls='none',
                         label="pFIRE")
                plt.title("Map {} component, R^2={:0.3}".format(dim, corr**2))
                plt.legend()
                plt.savefig(savepath)
                plt.close()
                image_entries.append(".. image:: {}"
                                     "".format(os.path.basename(savepath)))
        except IndexError:
            break

    print(tabulate(table_entries, headers="firstrow", tablefmt="grid"))

    table = tabulate(table_entries, headers="firstrow", tablefmt="rst")

    return (table, "\n".join(image_entries))
