#!/usr/bin/env python3

""" I/O routines for maps, masks and images
"""

from functools import partial
import numpy as np
import h5py
import sys
import skimage.io as skio

import flannel.io as fio

def __LINE__(): # FIXME remove
    return str(sys._getframe(1).f_lineno)


def load_image(imagepath):
    """ Attempt to load image with known image loaders
    """
    print("\n Load_image = {}".format(imagepath))
    print("imagepath" + imagepath)
    
    image = None
    imgloaders = [load_pfire_image, fio.load_image
                  ,partial(skio.imread, as_gray=True)
                  ]

    for loader in imgloaders:
        print(" trying loader="+ str(imgloaders.index(loader)))
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

    print("Loop over map loaders")
    print("mappath={}", mappath)
    
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
    print("\nload_pfire_image filepath={}".format(imagepath) )
    imagepath+=":/registered"  # FIXME data group should be read from xdmf file
    
    filename, group = [x.strip() for x in imagepath.split(':')]
    print("\nload_pfire_image filename={}".format(filename) )
    print("\nload_pfire_image group={}".format(group ))
    print("load_pfire_image() line="+__LINE__())
             
    # Workaround to bypass issue CBM-66 Xdmf lib bug. TO read metadata from XDMF
    if filename.endswith(".xdmf"):    
        filename += ".h5"
    
       
    print("\nload_pfire_image filename={}".format(filename) )
    print("load_pfire_image() line="+__LINE__())
    
    with h5py.File(filename, 'r') as fh:
        print(list(fh.keys()))
        imgdata = np.asarray(fh['registered'])

    return imgdata


def load_pfire_map(imagepath):
    """Load pFIRE map
    """
    imagepath+=":/map"  # FIXME
    
    print("\nload_pfire imagepath={}", imagepath)
    filename, group = [x.strip() for x in imagepath.split(':')]
    if filename.endswith(".xdmf"):
        filename += ".h5"

    print("\nload_pfire_map={}", load_pfire_map)
    
    print("\nload_pfire_filename={}", filename)
    
    imgdata = []
    with h5py.File(filename, 'r') as fh:
        for dxn in ['x', 'y', 'z']:
            print("\n dxn={}"+str(dxn))
            try:
                imgdata.append(np.asarray(fh["{}/{}".format(group, dxn)]))
            except KeyError:
                break

    return np.stack(imgdata)
