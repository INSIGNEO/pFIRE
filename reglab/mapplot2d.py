#!/usr/bin/env python3

import argparse

import numpy as np
import h5py

import matplotlib.figure as mf
import matplotlib.backends.backend_agg as mplbea

import scipy.interpolate as si

import skimage.io as skio

basesize = 8

def parse_args():
    parser = argparse.ArgumentParser(description="Create 2D png quiverplots "
                                                 "from hdf5 map data.")
    parser.add_argument("fixed", help="Path to fixed image")
    parser.add_argument("moved", help="Path to moved image")
    parser.add_argument("map", help="Path to map file")
    parser.add_argument("output", help="Path to save output png.")
    parser.add_argument("--group", default="map", help="Group name of map.")
    parser.add_argument("--flip", action="store_true", default=False,
                        help="Flip image vertically (origin top-left)")
    parser.add_argument("--invert-map", action="store_true", default=False,
                        help="Invert mapping before plotting")

    return parser.parse_args()


def reverse_mapping(nodes, mapping):
    assert(len(nodes)+1 == mapping.ndim)
    assert(tuple(len(x) for x in nodes) == mapping.shape[1:])

    nodegrid = np.stack(np.meshgrid(nodes[0], nodes[1], indexing='ij'))

    ofsnodes = np.moveaxis((nodegrid + mapping).reshape(nodegrid.shape[0], -1),
                           0, -1)
    rev_vecs = (-mapping).reshape(mapping.shape[0], -1)

    return np.stack(si.griddata(ofsnodes, rv, np.moveaxis(nodegrid, 0, -1),
        method="cubic") for rv in rev_vecs)


def get_map_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        shape = (2,) + h5f["{}/x".format(groupname)].shape
        data = np.empty(shape)
        for idx, tmpl in enumerate(("{}/x", "{}/y")):
            data[idx] = np.asarray(h5f[tmpl.format(groupname)])
    return data

def get_nodes_from_h5(mapfile, groupname):
    with h5py.File(mapfile, 'r') as h5f:
        data = []
        for tmpl in ("{}/nodes_x", "{}/nodes_y"):
            data.append(np.asarray(h5f[tmpl.format(groupname)]))
    return data

def main():
    args = parse_args()

    mapdata = get_map_from_h5(args.map, args.group)
    assert mapdata.ndim == 3

    nodes_x, nodes_y = get_nodes_from_h5(args.map, args.group)

    ns_x = np.diff(nodes_x)[0]/2
    ns_y = np.diff(nodes_y)[0]/2

    fixed = skio.imread(args.fixed, as_grey=True)
    moved = skio.imread(args.moved, as_grey=True)

    assert fixed.ndim == 2
    assert fixed.shape == moved.shape

    if args.invert_map:
        mapdata = reverse_mapping((nodes_x, nodes_y), mapdata)

    figsize = (basesize, basesize*(len(nodes_y)/len(nodes_x)))

    fig = mf.Figure(figsize=figsize)
    mplbea.FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    nnx, nny = np.meshgrid(nodes_x, nodes_y, indexing='ij')
    ax.imshow(fixed, origin='lower',
              extent=[0, fixed.shape[1], 0, fixed.shape[0]],
              cmap="Greys_r", alpha=1.0)
    ax.imshow(moved, origin='lower',
              extent=[0, moved.shape[1], 0, moved.shape[0]],
              cmap="Reds_r", alpha=0.5)
    ax.quiver(nnx, nny, mapdata[0], mapdata[1],
              angles='xy', scale_units='xy', scale=1)
    ax.set_xlim(nodes_x.min() - ns_x, nodes_x.max() + ns_x)
    ax.set_ylim(nodes_y.min() - ns_y, nodes_y.max() + ns_y)

    if args.flip:
        ax.invert_yaxis()

    fig.savefig(args.output, fmt="png")

    print("Saved {}".format(args.output))

if __name__ == "__main__":
    main()
