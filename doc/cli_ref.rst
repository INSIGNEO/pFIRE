===================
pFIRE Configuration
===================

Required Parameters
===================

These parameters must be provided, or pFIRE will abort with an error.

.. option:: fixed 

   Path to fixed image.

   May be absolute path, or relative to working directory.   

.. option:: moved 
   
   Path to moved image.

.. option:: nodespacing

   Desired final nodal spacing of displacement map.


Optional Parameters
===================

Default values exist for these parameters, but they may be overridden by the user.

.. option:: registered 

   *default: "registered.xmdf:/registered"*
   
   Path (and group name) to save registered image. 

   Output writer is chosed based on the file extension.  Current writers include hdf5, xdmf, and
   OpenImageIO (png and jpg - 2D only).

   Syntax is "path:group-path", where `group-path` sets the path and dataset name within the HDF5
   and dataset name within XDMF. `group-path` will be ignored by image writers which do not support
   it.

   Path may be an absolute path or relative to the working directory.  Group name may be any valid
   HDF5 group path.

.. option:: map 

   *default: "map.xdmf:/map"*
   
   Path (and group name) to save displacement map.

   Option syntax and behaviour is the same as for `registered` (see above).

   N.B image format writers do not support writing map data.

.. option:: verbose 

   *default: false*
   
   Boolean flag to enable verbose logging output.

.. option:: debug_frames 

   *default: false* 
   
   Output all intermediate image and map frames in the registration.
   
   Debug image frames will be saved to `%registered%-debug-%s-%i%` where `%s%` is the step
   number, `%i%` the iteration number, and `%registered%` the registered image output path. Map
   frames will be output following the same pattern using the map output path.
