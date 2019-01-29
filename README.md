# pFIRE
## The Parallel Framework for Image Registration

pFIRE performs elastic registration of 2- and 3-dimensional images by the method of Barber and
Hose&#91;[1](#note1)&#93;. It is implemented in C++ using the
[https://www.mcs.anl.gov/petsc/](PETSc) scientific toolkit to provide parallelised mathematical
routines.

Usage
-----

pFIRE is invoked from the console (via mpirun for parallel use), e.g:

  ```sh
  $ pfire config.cfg
  # or
  $ mpirun -np 8 pfire config.cfg
  ```

It takes a single argument - the name of a configuration file.  This is a deliberate design
decision to encourage documentation of the performed analysis.  At minimum this configuration file
must detail the path to the fixed image, path to the moved image, and required final spacing of the
map nodes. It also takes various optional arguments allowing customisation of the registration
process and outputs. pFIRE will automatically identify and open all supported file formats.  The
currently non-exhaustive list includes: `dicom`, ShIRT `.image` and `.mask` files, and the majority
of common 2-dimensional image formats via the OpenImageIO library.

A minimal usage example would be:

  ```ini
  fixed = fixed_image.dcm
  moved = "moved image.dcm"
  nodespacing = 10
  ```

This will register `fixed_image.dcm` to `moved image.dcm` with a nodespacing of 10 pixels, note
that if the filename contains spaces it must be enclosed in quotes.  In this case the output map
and registered image will be saved to the default hdf5+xdmf format with filenames
``registered.xdmf`` and ``map.xdmf``.

Dependencies
------------

pFIRE is an MPI application built on top the the [https://www.mcs.anl.gov/petsc/](PETSc)
distributed scientific toolkit. It is additionally dependent on the [https://www.boost.org/](Boost)
libraries for general utility routines and [https://www.hdfgroup.org/solutions/hdf5/](HDF5) input
and output support.  There are also optional dependencies which allow support for various
additional input image file formats.  At least one of these should be used depending on your
intended use case.

For all dependencies we recommend using the latest stable version.  Additionally we recommend that
[https://www.mcs.anl.gov/petsc/](PETSc) be configured to using single precision floating point
numbers.  This halves memory usage as double precision math provides no real benefit for image
registration.

**Required Dependencies**

   * [https://www.mcs.anl.gov/petsc/](PETSc) >= 3.10.0 (Recommend --with-precision=single)
   * [https://www.boost.org/](Boost) >= 1.58
   * [https://www.hdfgroup.org/solutions/hdf5/](HDF5) >= 1.10.0

*Optional Dependencies*

   * [https://dicom.offis.de/dcmtk.php.en](DCMTK) >= 3.6.3 (Support for DICOM image input)
   * [http://www.openimageio.org/](OpenImageIO) >= 1.8.13 (General purpose image format support e.g .png .tiff and image stack support)

We recommend installing dependencies using your system package manager (e.g synaptic, apt, yum), or
on HPC the use of SPACK may be appropriate.

Installing
----------
pFIRE is installed using cmake in an out of tree build:

```sh
user@machine $ tar xvf pfire-latest.tar.gz
user@machine $ cd pfire
user@machine $ mkdir build
user@machine $ cd build
user@machine $ cmake ..
user@machine $ make
```

Links
-----
<a name="note1">[1]</a>: DC Barber and DR Hose 2005 (https://doi.org/10.1080/03091900412331289889),
    DC Barber et al. 2007 (https://doi.org/10.1016/j.media.2007.06.011)

<a name="note2">[2]</a>: https://www.mcs.anl.gov/petsc/
