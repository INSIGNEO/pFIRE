#!/usr/bin/env python3

""" I/O routines for maps, masks and images
"""

from functools import partial
import numpy as np
import h5py

import skimage.io as skio

import flannel.io as fio


def load_image(imagepath):
    """ Attempt to load image with known image loaders
    """
    image = None
    imgloaders = [load_pfire_image, fio.load_image,
                  partial(skio.imread, as_gray=True)]
    for loader in imgloaders:
        try:
            image = loader(imagepath)
        except (ValueError, OSError):
            pass
        if image is not None:
            break

    if image is None:
        raise RuntimeError("Failed to load image \"{}\"".format(imagepath))

    image = image / image.max()

    return image


def load_map(mappath):
    """ Attempt to load map with known loaders
    """
    data = None
    shirtloader = lambda path: fio.load_map(path)[0][0:3]
    maploaders = [load_pfire_map, shirtloader]

    for loader in maploaders:
        try:
            data = loader(mappath)
        except (ValueError, OSError):
            pass
        if data is not None:
            break

    if data is None:
        raise RuntimeError("Failed to load map \"{}\"".format(mappath))

    return data

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
