#!/usr/bin/env python3

import numpy as np

def save_image(im, fname):
    sizearr = np.ones(4, dtype=np.int32)
    sizearr[0:len(im.shape)] = im.shape
    im = im.astype(np.float32).flatten(order='F')
    with open(fname, 'wb') as fh:
        fh.write(sizearr.tobytes())
        fh.write(im.tobytes())


def load_image(fname):
  with open(fname, 'rb') as fh:
    shape = np.frombuffer(fh.read(16), np.int32)
    imsize = np.product(shape) * np.float32(0).itemsize
    im = np.frombuffer(fh.read(imsize), np.float32)
  im = np.squeeze(im.reshape(shape, order='F'))
  return im


def save_mask(im, fname):
  sizearr = np.ones(3,dtype=np.int32)
  sizearr[0:len(im.shape)] = im.shape
  im = im.astype(np.int16).flatten(order='F')
  with open(fname, 'wb') as fh:
    fh.write(sizearr.tobytes())
    fh.write(im.tobytes())
    

def load_mask(fname):
  with open(fname, 'rb') as fh:
    shape = np.frombuffer(fh.read(12), np.int32)
    imsize = np.product(shape) * np.float32(0).itemsize
    im = np.frombuffer(fh.read(imsize), np.float32)
  im = np.squeeze(im.reshape(shape, order='F'))
  return im


def save_map(im, grid, fname, imgdims=None):
  """
  Takes map as (N,X,..) matrix and creates shirt mapfile using
  parameter grid to map node locations
  """
  if grid.shape != im.shape[1:]:
    raise ValueError("Dimensionality mismatch imgmap vs grid")

  if imgdims is not None:
    if len(imgdims) != len(im.shape):
      raise ValueError("Dimensionality mismatch imgdims vs imgmap.shape")

  shape3d = np.ones(3)
  shape3d[0:len(grid.shape)] = grid.shape
  gridmesh = np.mgrid[0:shape3d[0]
                     ,0:shape3d[1]
                     ,0:shape3d[2]]

  len1d = np.product(shape3d)

  mapheader = np.zeros(9, np.int32)
  mapheader[0:3] = shape3d #map dimensions
  mapheader[3] = len1d #total data size
  mapheader[4] = gridmesh[0][1,0,0] - gridmesh[0][0,0,0] #node spacing
  mapheader[5] = grid.shape # dimensionality
  mapheader[6:9] = imgdims # image dimensions

  mapdata = np.zeros((7,len1d), np.int32)

  for i in range(3):
    mapdata[i] = gridmesh[i].reshape(-1)

  if im.shape[-1] == 3:
    mapdata[3] = im[0]
    mapdata[6] = im[-2]
  if im.shape[-1] == 4:
    mapdata[3] = im[0]
    mapdata[4] = im[1]
    mapdata[6] = im[-2]
  if im.shape[-1] == 5:
    mapdata[3] = im[0]
    mapdata[4] = im[1]
    mapdata[5] = im[2]
    mapdata[6] = im[-2]

  with open(fname, 'wb') as fh:
    fh.write(mapheader.tobytes()) 
    fh.write(mapdata.tobytes()) 


def load_map(fname):
  """
  Loads shirt map file returning map as (N,X,..) matrix and grid as a list of N
  vectors
  """
  with open(fname, 'rb') as fh:
    mapheader = np.frombuffer(fh.read(9*np.int32(4).itemsize), np.int32)
    mapsize = tuple(mapheader[0:3])
    mapdatasize = int(np.float64(0).itemsize * np.product(mapsize))
    mapdata = []
    for i in range(7):
      databuf = np.frombuffer(fh.read(mapdatasize), np.float64)
      mapdata.append(databuf.reshape(mapsize, order='F'))

  griddata = [mapdata[0][:,0,0], mapdata[1][0,:,0], mapdata[2][0,0,:]]
  
  mapdata = np.stack(mapdata[3:])

  return(mapdata, griddata)
