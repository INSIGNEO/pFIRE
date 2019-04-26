#!/usr/bin/env python3

import numpy as np
import h5py

def load_pfire_image(imagepath):
    """Load pFIRE image
    """
    filename, group = [x.strip() for x in imagepath.split(':')]
    if filename.endswith(".xdmf"):
        filename += ".h5"

    with h5py.File(filename, 'r') as fh:
        imgdata = np.asarray(fh[group])

    return imgdata


def load_pfire_map(imagepath):
    """Load pFIRE map
    """
    filename, group = [x.strip() for x in imagepath.split(':')]
    if filename.endswith(".xdmf"):
        filename += ".h5"

    imgdata = []
    with h5py.File(filename, 'r') as fh:
        for dxn in ['x', 'y', 'z']:
            try:
                imgdata.append(np.asarray(fh["{}/{}".format(group, dxn)]))
            except KeyError:
                break

    return np.stack(imgdata)


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
