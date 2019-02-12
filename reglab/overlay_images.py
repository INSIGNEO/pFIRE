#!/usr/bin/env python3

import argparse

import numpy as np

import matplotlib.figure as mf
import matplotlib.backends.backend_agg as mplbea

import skimage.io as skio

dpi=96

def parse_args():
    parser = argparse.ArgumentParser(description="Overlay images with "
                                     "transparency")
    parser.add_argument("fixed", help="Path to fixed image")
    parser.add_argument("moved", help="Path to moved image")
    parser.add_argument("output", help="Path to save output png.")
    parser.add_argument("--flip", action="store_true", default=False,
                        help="Flip image vertically (origin top-left)")

    return parser.parse_args()

def main():
    args = parse_args()

    fixed = skio.imread(args.fixed, as_grey=True)
    moved = skio.imread(args.moved, as_grey=True)

    assert fixed.ndim == 2
    assert fixed.shape == moved.shape

    figsize = (fixed.shape[1]/dpi, fixed.shape[0]/dpi)

    fig = mf.Figure(figsize=figsize)
    mplbea.FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.imshow(fixed, origin='lower',
              extent=[0, fixed.shape[1], 0, fixed.shape[0]],
              cmap="Greys_r", alpha=1.0)
    ax.imshow(moved, origin='lower',
              extent=[0, moved.shape[1], 0, moved.shape[0]],
              cmap="Reds_r", alpha=0.5)

    if args.flip:
        ax.invert_yaxis()

    fig.savefig(args.output, fmt="png", dpi=dpi)

    print("Saved {}".format(args.output))

if __name__ == "__main__":
    main()
