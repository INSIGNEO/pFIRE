#!/usr/bin/env python3

import os
import re

import requests
from tqdm import tqdm

import numpy as np
import skimage.io as skio
import skimage.transform as skt
import skimage.exposure as ske

base_url = "http://ptak.felk.cvut.cz/Medical/dataset_CIMA/lung-lesion_1/scale-10pc"

image_list = ["29-041-Izd2-w35-CD31-3-les1",
              "29-041-Izd2-w35-Cc10-5-les1",
              "29-041-Izd2-w35-He-les1",
              "29-041-Izd2-w35-Ki67-7-les1",
              "29-041-Izd2-w35-proSPC-4-les1"]

temp_tmpl = "{}-temp"
img_tmpl = "{}.jpg"
csv_tmpl = "{}.csv"

download_auth = ('dataset-guest', 'icip2018')

prere = re.compile("29-041-Izd2-w35-")
postre = re.compile("-.*")

target_shape=(850,700)

def main():
    for image in image_list:
        fetch_base = "{}/{}".format(base_url, image)
        output_base = prere.sub("lung_lesion_", image)
        output_base = postre.sub("", output_base).lower()
        output_temp = temp_tmpl.format(output_base)
        imgtmp, csvtmp = get_image_pair(fetch_base, output_temp, download_auth)

        print("Resizing {} to {}x{}".format(imgtmp, *target_shape))
        imgdata = skio.imread(imgtmp, as_grey=True)
        imgdata = ske.equalize_adapthist(imgdata)
        full_shape = imgdata.shape
        imgdata = skt.resize(imgdata, target_shape, mode="edge")
        skio.imsave(img_tmpl.format(output_base), imgdata)
        csvdata = np.loadtxt(csvtmp, skiprows=1, dtype=np.int32, delimiter=",")
        csvdata[:, 1] = csvdata[:, 1]/(full_shape[1]/target_shape[1])
        csvdata[:, 2] = csvdata[:, 2]/(full_shape[0]/target_shape[0])
        csvdata = csvdata.astype(np.int32)
        np.savetxt(csv_tmpl.format(output_base), csvdata, header=",X,Y",
                   delimiter=",", fmt="%i")


def get_image_pair(fetch_base, output_base, auth):
    paths = []
    for tmpl in (img_tmpl, csv_tmpl):
        pth = tmpl.format(output_base)
        pretty_get(tmpl.format(fetch_base), pth, auth)
        paths.append(pth)
    return paths

def pretty_get(url, output_path, authdata=None):
    response = requests.get(url, stream=True, auth=authdata)
    block_total = int(response.headers.get('content-length', 0))
    try:
        if(os.path.getsize(output_path) == block_total):
            print("File {} already fetched".format(output_path))
            return
    except:
        pass
    with open(output_path, "wb") as fh:
        block_size = 1024*1024
        with tqdm(total=block_total, unit='B', unit_scale=True, leave=False) as pb:
            for data in response.iter_content(chunk_size=block_size):
                pb.update(block_size)
                fh.write(data)

if __name__ == "__main__":
    main()
