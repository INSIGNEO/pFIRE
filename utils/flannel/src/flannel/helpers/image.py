#!/usr/bin/env python3

import sys
import os.path
import skimage.io as skio
import numpy as np
import pydicom

from flannel.io import load_image, save_image, save_mask

def image_to_shirt():
    if len(sys.argv) != 2:
        print("Usage: {} /path/to/imagefile.ext".format(sys.argv[0]))
        print("\tOutput: imagefile.image")
        return(0)

    input_path = sys.argv[1]
    input_dir, input_file = os.path.split(input_path)
    output_file = os.path.splitext(input_file)[0] + '.image'
    output_path = output_file

    try:
        image = skio.imread(input_path, as_grey=True).astype(np.float32)
    except Exception as err:
        try:
            dcm = pydicom.dcmread(input_path)
            image = dcm.pixel_array
        except Exception as err:
            print("Error: Failed to open {}".format(input_path))
            print(err)
            return(-1)

    save_image(image, output_path)

def image_to_mask():
    if len(sys.argv) != 2:
        print("Usage: {} /path/to/imagefile.ext".format(sys.argv[0]))
        print("\tOutput: imagefile.mask")
        return(0)

    input_path = sys.argv[1]
    input_dir, input_file = os.path.split(input_path)
    output_file = os.path.splitext(input_file)[0] + '.mask'
    output_path = output_file

    try:
        image = skio.imread(input_path, as_grey=True)
        mask = np.ones(image.shape, dtype=np.int16)
        save_mask(mask, output_path)
        
    except Exception as err:
        try:
            dcm = pydicom.dcmread(input_path)
            image = dcm.pixel_array
            mask = np.ones(image.shape, dtype=np.int16)
            save_mask(mask, output_path)
        except Exception as err:
            print("Error: Failed to open {}".format(input_path))
            print(err)
            return(-1)

        #save_mask(mask, output_path)


def shirt_to_image():
    if len(sys.argv) != 2:
        print("Usage: {} /path/to/imgname.image".format(sys.argv[0]))
        print("\tOutput: imagename.png")
        return(0)

    input_path = sys.argv[1]
    input_dir, input_file = os.path.split(input_path)
    output_file = os.path.splitext(input_file)[0] + '.png'
    output_path = output_file

    try:
        image = load_image(input_path)
    except Exception as err:
        print("Error: Failed to open {}".format(input_path))
        print(err)
        return(-1)

    print(image.shape)
    print(image.min())
    print(image.max())

    image = image - image.min()
    image = image / image.max()

    skio.imsave(output_path, image)
