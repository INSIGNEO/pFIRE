#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import scipy.interpolate as si
import h5py

import skimage.io as skio

import matplotlib.figure as mf
import matplotlib.backends.backend_agg as mplbea

cols1 = ('Greys_r', 'black')
cols2 = ('Reds_r', 'red')

dpi = 300

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(-1)

def parse_args():
    parser = MyParser(description="Plot annotated histology slides with map "
                                  "displacements")
    parser.add_argument('fixed_img', help="Fixed image")
    parser.add_argument('fixed_annot', help="Fixed image annotation csv")
    parser.add_argument('moved_img', help="Moved image")
    parser.add_argument('moved_annot', help="Moved image annotation csv")
    parser.add_argument('map', help="pFIRE Map Output")
    parser.add_argument('output', help="Output image")
    parser.add_argument('--map_group', default="map", help="Map hdf5 group")

    args = parser.parse_args()

    return args

def get_map_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        npshape = (2,) + h5f["{}/x".format(groupname)].shape
        data = np.zeros(npshape)
        for idx, tmpl in enumerate(("{}/x", "{}/y")):
            data[idx,:,:] = h5f[tmpl.format(groupname)]
    return data


def get_nodes_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        data = []
        for tmpl in ("{}/nodes_x", "{}/nodes_y"):
            data.append(np.asarray(h5f[tmpl.format(groupname)]))
    return data

    mapdata = get_map_from_h5(args.map, args.map)

def main():

    args = parse_args()
    
    fixed_img = skio.imread(args.fixed_img, as_gray=True)
    fixed_annot = np.loadtxt(args.fixed_annot, delimiter=',', skiprows=1, dtype=np.int32)

    moved_img = skio.imread(args.moved_img, as_gray=True)
    moved_annot = np.loadtxt(args.moved_annot, delimiter=',', skiprows=1, dtype=np.int32)

    # load map data
    mapdata = get_map_from_h5(args.map, args.map_group)
    nodes_x, nodes_y = get_nodes_from_h5(args.map, args.map_group)

    # update csv points from interpolated map
    interp_x = si.RegularGridInterpolator((nodes_x, nodes_y), mapdata[0])
    interp_y = si.RegularGridInterpolator((nodes_x, nodes_y), mapdata[1])

    offsets_x = interp_x(fixed_annot[:,1:])
    offsets_y = interp_y(fixed_annot[:,1:])

    # Set figsize for 1:1 pixels
    figsize = (fixed_img.shape[1]/dpi, fixed_img.shape[0]/dpi)

    fig = mf.Figure(figsize=figsize)
    canvas = mplbea.FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.imshow(fixed_img, cmap=cols1[0], aspect='auto')
    ax.imshow(moved_img, cmap=cols2[0], alpha=0.5, aspect='auto')
    ax.plot(fixed_annot[:, 1], fixed_annot[:, 2], color=cols1[1], marker='.',
            linestyle="none", markersize="1")
    ax.plot(moved_annot[:, 1], moved_annot[:, 2], color=cols2[1], marker='.',
            linestyle="none", markersize="1")

    ax.quiver(fixed_annot[:, 1], fixed_annot[:, 2], offsets_x, offsets_y,
              angles='xy', scale_units='xy', scale=1)
    fig.savefig(args.output, dpi=dpi)


if __name__ == "__main__":
    main()
