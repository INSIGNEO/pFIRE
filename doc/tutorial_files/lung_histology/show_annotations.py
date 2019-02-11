#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np

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
    parser = MyParser(description="Plot and optionally compare annotated histology images")
    parser.add_argument('output', help="Output image")
    parser.add_argument('img1', help="First image")
    parser.add_argument('annot1', help="First annotation csv")
    parser.add_argument('other_img', nargs=argparse.REMAINDER, help="Second image and csv")

    args = parser.parse_args()

    if len(args.other_img) == 0:
        pass
    elif len(args.other_img) == 2:
        vars(args)['img2'] = args.other_img[0]
        vars(args)['annot2'] = args.other_img[1]
    else:
        parser.error("Can specify exactly one extra image+annotations")

    return args


def main():

    args = parse_args()
    
    img1 = skio.imread(args.img1, as_gray=True)
    annot1 = np.loadtxt(args.annot1, delimiter=',', skiprows=1, dtype=np.int32)

    second_img = True
    try: 
        img2 = skio.imread(args.img2, as_gray=True)
        annot2 = np.loadtxt(args.annot2, delimiter=',', skiprows=1, dtype=np.int32)
    except:
        second_img = False

    figsize = (img1.shape[1]/dpi, img1.shape[0]/dpi)

    fig = mf.Figure(figsize=figsize)
    canvas = mplbea.FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.imshow(img1, cmap=cols1[0], aspect='auto')
    if second_img:
        ax.imshow(img2, cmap=cols2[0], alpha=0.5, aspect='auto')
    ax.plot(annot1[:, 1], annot1[:, 2], color=cols1[1], marker='.',
            linestyle="none", markersize="1")
    if second_img:
        ax.plot(annot2[:, 1], annot2[:, 2], color=cols2[1], marker='.',
                linestyle="none", markersize="1")

    fig.savefig(args.output, dpi=dpi)


if __name__ == "__main__":
    main()
