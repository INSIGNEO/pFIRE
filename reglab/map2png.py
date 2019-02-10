#!/usr/bin/env python3

import argparse

import numpy as np
import h5py

import matplotlib.figure as mf
import matplotlib.backends.backend_agg as mplbea

import skimage.io as skio

basesize = 8

def parse_args():
    parser = argparse.ArgumentParser(description="Create 2D png quiverplots "
                                                 "from hdf5 map data.")
    parser.add_argument("fixed", help="Path to fixed image")
    parser.add_argument("moved", help="Path to moved image")
    parser.add_argument("map", help="Path to map file")
    parser.add_argument("output", help="Path to save quiverplot png.")
    parser.add_argument("--group", default="map", help="Group name of map.")
    parser.add_argument("--flip", action="store_true", default=False)

    return parser.parse_args()


def get_map_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        data = []
        for idx, tmpl in enumerate(("{}/x", "{}/y")):
            data.append(np.asarray(h5f[tmpl.format(groupname)]))

    return data

def get_nodes_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        data = []
        for idx, tmpl in enumerate(("{}/nodes_x", "{}/nodes_y")):
            data.append(np.asarray(h5f[tmpl.format(groupname)]))

    return data


def main():
    args = parse_args()
    
    nodes_x, nodes_y = get_nodes_from_h5(args.map, args.group) 
    data_x, data_y = get_map_from_h5(args.map, args.group) 
    
    fixed = skio.imread(args.fixed, as_gray=True)
    moved = skio.imread(args.moved, as_gray=True)

    figsize = (basesize, basesize*(data_x.shape[1]/data_x.shape[0]))

    fig = mf.Figure(figsize=figsize)
    mplbea.FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    nnx, nny = np.meshgrid(nodes_x, nodes_y, indexing='ij')
    ax.imshow(fixed, origin='lower',
              extent=[0, fixed.shape[0], 0, fixed.shape[1]],
              cmap="Reds_r", alpha=0.5)
    ax.imshow(moved, origin='lower',
              extent=[0, moved.shape[0], 0, moved.shape[1]],
              cmap="Reds_r", alpha=0.5)
    ax.quiver(nnx, nny, data_x, data_y)
    print(nnx.min(), nnx.max())
    print(fixed.shape)
    ax.set_xlim(nodes_x.min(), nodes_x.max())
    print(nodes_x.min(), nodes_x.max())
    ax.set_ylim(nodes_y.min(), nodes_y.max())
    print(nodes_y.min(), nodes_y.max())

    if(args.flip):
        ax.invert_yaxis()

    fig.savefig(args.output)

if __name__ == "__main__":
    main()
