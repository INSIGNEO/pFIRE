#!/usr/bin/env python3

import numpy as np

from flannel import io as fio

image_shape = (140, 100, 60)
block_shape = (50, 40, 40)
block1_fixed_offset = (10, 10, 10)
block2_fixed_offset = (80, 10, 10)
block1_moved_offset = (10, 50, 10)
block2_moved_offset = (80, 50, 10)
mask_shape = (70, 100, 60)
mask_offset = (0, 0, 0)

def slice_from_offsets(shape, offset):
    return tuple(slice(x, x+y) for x, y in zip(offset, shape))


def main():
    """Create 3D masking test data
    """

    print("Saving fixed.image")
    fixed = np.zeros(image_shape)
    fixed[slice_from_offsets(block_shape, block1_fixed_offset)] = 1
    fixed[slice_from_offsets(block_shape, block2_fixed_offset)] = 1
    fio.save_image(fixed, "fixed.image")

    print("Saving moved.image")
    moved = np.zeros(image_shape)
    moved[slice_from_offsets(block_shape, block1_moved_offset)] = 1
    moved[slice_from_offsets(block_shape, block2_moved_offset)] = 1
    fio.save_image(moved, "moved.image")

    print("Saving mask.mask")
    mask = np.zeros(image_shape)
    mask[slice_from_offsets(mask_shape, mask_offset)] = 1
    fio.save_mask(mask, "mask.mask")

if __name__ == "__main__":
    main()
