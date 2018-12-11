=============
Running pFIRE
=============

Interaction with pFIRE on the command line is achieved through the `pfire` executable. There are
various subcommands available to allow interaction with images, maps, masks and various
registration options.  All the available options are documented in the :doc:`command
reference<cli_ref>`.

Elastic Registration
--------------------

pFIRE is invoked from the console:

.. runblock:: console

  $ pfire --help

It takes a single argument - the name of a configuration file.  This is a deliberate design
decision to encourage documentation of the performed analysis.  At minimum this configuration file
must detail the path to the fixed image, path to the moved image, and required final spacing of the
map nodes. It also takes various optional arguments allowing customisation of the registration
process and outputs. pFIRE will automatically identify and open all supported file formats.  The
currently non-exhaustive list includes: `dicom`, ShIRT `.image` and `.mask` files, and the majority
of common 2-dimensional image formats via the OpenImageIO library.

A minimal usage example would be:

.. code-block:: ini

  fixed = fixed_image.dcm
  moved = "moved image.dcm"
  nodespacing = 10

This will register `fixed_image.dcm` to `moved image.dcm` with a nodespacing of 10 pixels, note
that if the filename contains spaces it must be enclosed in quotes.  In this case the output map
and registered image will be saved to the default hdf5 format with filename ``registration.h5``.


ShIRT Compatibility
-------------------

pFIRE is designed to be compatible with ShIRT both algorithmically and practially.  Therefore we
provide a way to use pFIRE as a drop in replacement for ``ShIRT Register`` in existing workflows.
If the executable is invoked with the name `shirt` or `pfire-shirt` (either by renaming or using a
symlink) then pFIRE will interpret command line input using the ShIRT syntax.  In this mode the
output formats are fixed as ShIRT format for use with existing analysis workflows.

.. code-block:: shell

  # symlink example
  $ ln -s ./pfire ./pfire-shirt
  $ pfire-shirt register fixed fixed.image moved moved.image nodespacing 10

  # copy example
  $ cp -s ./pfire ./shirt
  $ shirt register fixed fixed.image moved moved.image mask mask.mask nodespacing 10
