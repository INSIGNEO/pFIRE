#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import h5py
import scipy.interpolate as si


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(-1)

def parse_args():
    parser = MyParser(description="Update annotations with map data")
    parser.add_argument('annot', help="Input annotation csv")
    parser.add_argument('map', help="Input annotation csv")
    parser.add_argument('nodespacing', type=float, help="node spacing")
    parser.add_argument('image_shape', type=int, nargs=2, help="image spacing (space separated e.g 200 300)")
    parser.add_argument('output', help="Output annotation csv")
    parser.add_argument('--map_name', default="map", help="HDF5 group name")

    args = parser.parse_args()

    return args

def get_map_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        npshape = (2,) + h5f["{}/x".format(groupname)].shape
        data = np.zeros(npshape)
        for idx, tmpl in enumerate(("{}/x", "{}/y")):
            data[idx,:,:] = h5f[tmpl.format(groupname)]

    return data

def main():

    args = parse_args()

    # load csv data
    csvdata = np.loadtxt(args.annot, delimiter=',', skiprows=1,
                         dtype=np.float32)

    # load map data
    mapdata = get_map_from_h5(args.map, args.map_name)
    mapdata = np.flip(mapdata, 1)

    # calculate node locations:
    nodelen_x =args.nodespacing * (mapdata.shape[1] - 1) 
    nodelen_y =args.nodespacing * (mapdata.shape[2] - 1) 
    offset_x = 0.5*(nodelen_x - args.image_shape[0])
    offset_y = 0.5*(nodelen_y - args.image_shape[1])

    nodes_x = np.linspace(offset_x, offset_x + nodelen_x, mapdata.shape[1])
    nodes_y = np.linspace(offset_y, offset_y + nodelen_y, mapdata.shape[2])

    # update csv points from interpolated map
    interp_x = si.RegularGridInterpolator((nodes_x, nodes_y), mapdata[0])
    interp_y = si.RegularGridInterpolator((nodes_x, nodes_y), mapdata[1])

    offsets_x = interp_x(csvdata[:,1:])
    offsets_y = interp_y(csvdata[:,1:])

    csvdata[:,1] -= offsets_x
    csvdata[:,2] -= offsets_y

    import IPython; IPython.embed()

    # save csv data
    np.savetxt(args.output, csvdata, header=",X,Y", delimiter=",", fmt="%i")

if __name__ == "__main__":
    main()
