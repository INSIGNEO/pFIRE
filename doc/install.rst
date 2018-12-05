================
Installing pFIRE
================

Currently pFIRE must be compiled from source.  This is available from the github repository at
https://github.com/INSIGNEO/pfire_petsc.git

Dependencies
------------

pFIRE is an MPI application built on top the the PETSc_ distributed scientific toolkit. It is
additionally dependent on the Boost_ libraries for general utility routines and HDF5_ input and
output support.  There are also optional dependencies which allow support for various additional
input image file formats.  At least one of these should be used depending on your intended use
case.

For all dependencies we recommend using the latest stable version.  Additionally we recommend that
PETSc be configured to using single precision floating point numbers.  This halves memory usage as
double precision math provides no real benefit for image registration.

**Required Dependencies**

   * PETSC_ >= 3.10.0 (Recommend --with-precision=single)
   * Boost_ >= 1.58
   * HDF5_ >= 1.10.0

*Optional Dependencies*

   * DCMTK_ >= 3.6.3 (Support for DICOM image input)
   * OpenImageIO_ >= 1.8.13 (General purpose image format support e.g .png .tiff and image stack support)

We recommend installing dependencies using your system package manager (e.g synaptic, apt, yum), or
on HPC the use of SPACK_ may be appropriate.

.. _PETSc: https://www.mcs.anl.gov/petsc/
.. _Boost: https://www.boost.org/
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _DCMTK: https://dicom.offis.de/dcmtk.php.en
.. _OpenImageIO: http://www.openimageio.org/
.. _SPACK: https://spack.io


Building pFIRE
--------------

pFIRE is installed using cmake, so after checking out the code and ensuring all dependencies are
installed first call `cmake` followed by `make`

.. code-block:: shell

   git clone https://github.com/INSIGNEO/pfire_petsc.git
   cd pfire
   mkdir build
   cd build
   cmake ..
   make

