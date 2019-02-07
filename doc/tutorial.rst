==============
pFIRE Tutorial
==============

This tutorial is designed to provide a hands-on walkthrough of the capabilities of pFIRE.  We will
begin with an example that uses the default configuration of pFIRE.  After this we will look at
some of the more advanced options that allow us to control both how pFIRE performs its
registration, and the output formats that can be written.

All the tutorial files can be found on the `pFIRE Github`_, or may be downloaded in a :download:`zip archive
</_static/pfire_tutorial_files.zip>`.

.. _pfire github: https://github.com/INSIGNEO/pFIRE/tree/master/doc/tutorial_files

In order to perform registration, pFIRE requires a pair of images: a reference image, and a second
image that is deformed until it is the same as the reference image. These are referred to as the
**fixed** and **moved** images respectively.  These are set in the configuration file using the
keys ``fixed`` and ``moved``.

When performing the registration the applied deformation field is initially very coarse grained,
with just 3 grid nodes per dimension.  This is repeatedly refined to increase the number of nodes
in the registration until the required resolution is reached.  This is expressed in the
configuration file using the parameter ``nodespacing``, which gives the spacing between nodes in
image pixels.

A Minimal Configuration File
============================

pFIRE is run from the command line, with all the configuration options set in a configuration file.
This file uses ``.ini`` syntax, with ``key = value`` pairs.

The minimum information pFIRE requires to perform a registration is a fixed image, a moved image,
and a target nodespacing. If the parameters ``fixed``, ``moved`` and ``nodespacing`` are not set in
the configuration file, pFIRE will abort with an error. 

A minimal configuration file will look like this:

.. code-block:: ini

   fixed = path/to/fixed.img
   moved = path/to/moved.img
   nodespacing = 10

Example: Cartoon Faces
----------------------

This example is located in the ``faces_1`` directory, here you will find the files ``happy.png``,
``sad.png`` and ``sad2happy_default.conf``.  If you open ``sad2happy_default.conf`` in an editor
you will find the following:

.. literalinclude:: /tutorial_files/faces_1/sad2happy_default.conf
   :language: ini

This instructs pFIRE to register the image in ``sad.png`` to the image in ``happy.png`` using a
final nodespacing of 10.  The ``registered`` parameter instructs pFIRE to store the registered
image in the file ``sad2happy_default.png``.

+--------------------+--------------------+
| |happy|            | |sad|              |
+--------------------+--------------------+
| ``happy.png``      | ``sad.png``        |
+--------------------+--------------------+

This example can be run with:

.. code-block:: sh

   user@machine $ pfire sad2happy_default.conf

or using MPI with X tasks:

.. code-block:: sh

   user@machine % mpirun -np X sad2happy_default.conf

This will produce a pair of output files: ``sad2happy_default.png`` and ``map.xdmf``.  We will look
at the map output later in the tutorial. For now, viewing ``sad2happy_default.png`` shows the
results of the registration, with ``happy.png`` for comparison:

+----------------------------+----------------------------+
| |sad2happy|                | |happy|                    |
+----------------------------+----------------------------+
| ``sad2happy_default.png``  | ``happy.png``              |
+----------------------------+----------------------------+

If the registration was perfect ``sad2happy_default.png`` would be identical to ``happy.png``, however, we
can see that there a still differences: the mouth is somewhat distorted, and there are changes to
the eyes in the registered image even though they are identical in both the source images. This is
a result of the smoothing behaviour of pFIRE, which tries to ensure that the calculated
displacement is globally smooth. This is typically the desired behaviour, but in certain situations
the default smoothing behaviour is non-optimal for the problem at hand. To help with this, pFIRE
offers several configuration parameters to help control the smoothing behaviour.

.. |happy| image:: /tutorial_files/faces_1/happy.png

.. |sad| image:: /tutorial_files/faces_1/sad.png

.. |sad2happy| image:: /tutorial_files/faces_1/sad2happy_default.png


Customizing the Registration
============================

pFIRE has two key parameters that allow the user to optimize the registration, these are the value
of the smoothing parameter :math:`\lambda`, and whether or not the memory term is used in the
registration.  Both of these parameters have an effect on how smooth the displacement field is.

The Smoothing Parameter :math:`\lambda`
---------------------------------------

The registration equation [ref_] includes a smoothing constraint through the laplacian matrix,
which imposes a requirement for smoothness on the solution to the equation, and therefore on the
displacement field.  The value of the parameter :math:`\lambda` determines the relative strength of
the smoothing constraint relative to the registration constraint.

By default the value of :math:`\lambda` is automatically calculated such that the condition number
of the registration matrix :math:`\mathbf{T}^t\mathbf{T} + \lambda\mathbf{L}^t\mathbf{L}` is
minimized, however, in certain situations it may be useful to customize the smoothing behaviour.

The value and behaviour of :math:`\lambda` can be controlled by two configuration options:
``lambda_mult``, which allows the user to provide a scaling multiplier for the automatically
determined value of :math:`\lambda`; and ``lambda``, which allows the user to specify a fixed value
of :math:`\lambda`.

The ``lambda`` parameter
""""""""""""""""""""""""

The default value of the ``lambda`` parameter is ``auto``, if a value of ``lambda`` is not
specified in the configuration file then it defaults to this value.  If ``lambda = auto`` then the
value of lambda is calculated to minimize the condition number of the matrix 
:math:`\mathbf{T}^t\mathbf{T} + \lambda\mathbf{L}^t\mathbf{L}` as this maximises the robustness of
the algorithm.

An example configuration file for running pFIRE with a fixed value of :math:`\lambda` is provided
in the file ``faces_1/sad2happy_fixed_lambda.conf``.

.. literalinclude:: /tutorial_files/faces_1/sad2happy_fixed_lambda.conf
   :language: ini

This will register the same images as above, but with a large fixed ``lambda = 10``.

+---------------------------------+---------------------------------+
| |sad2happy_fixed_lambda|        | |happy|                         |
+---------------------------------+---------------------------------+
| ``sad2happy_fixed_lambda.png``  | ``happy.png``                   |
+---------------------------------+---------------------------------+

Setting ``lambda = 10`` in this way causes the algorithm to always give the smoothing 10 times the
weight of the registration when solving for the displacement field, and so in this case oversmooths
the image resulting in a poor registration. Additionally, fixing lambda in this way potentially
causes the registration matrix to be ill-conditioned and so can cause the solver to be unstable.
This functionality is included as it is sometimes appropriate for advanced users, however its use
in not generally recommended.  If the results are over- or under-smoothed, the ``lambda_mult``
parameter can be used to make relative adjustments to the smoothing.

.. |sad2happy_fixed_lambda| image:: /tutorial_files/faces_1/sad2happy_fixed_lambda.png

The ``lambda_mult`` parameter
"""""""""""""""""""""""""""""

The ``lambda_mult`` parameter is intended to be used in conjunction with ``lambda = auto`` to
allow the amount of smoothing to be adjusted relative to the calculated optimum.  This allows the
user some control of the amount of smoothing whilst retaining reasonable stability of the
registration algorithm.

In the first example above we saw that the mouth was not precisely registered because the image is
smoothed too strongly, therefore we might wish to try reducing the smoothing by setting the
parameter ``lambda_mult = 0.5``.

An example configuration file for this is provided in ``faces_1/sad2happy_relative_lambda.conf``:

.. literalinclude:: /tutorial_files/faces_1/sad2happy_relative_lambda.conf
   :language: ini

This registers the same images once again, but reducing the smoothing by a factor of 0.5:

+-----------------------------------+-----------------------------------+
| |sad2happy_relative_lambda|       | |happy|                           |
+-----------------------------------+-----------------------------------+
| ``sad2happy_relative_lambda.png`` | ``happy.png``                     |
+-----------------------------------+-----------------------------------+

.. |sad2happy_relative_lambda| image:: /tutorial_files/faces_1/sad2happy_relative_lambda.png

This results in a registration that is very similar to the first, the registration of the details
of the mouth is improved, but there are still issues with the rest of the image. This is a problem
with the way the image is smoothed, and cannot be solved by simply changing the smoothing strength,
instead we must change the fundamental smoothing behaviour.

The Memory Term
---------------

The registration algorithm implemented by pFIRE is an iterative one which applies the same equation
repeatedly to incrementally improve the displacement field until the registration is complete.  The
equation includes an optional extra smoothing term which depends on the total displacement field
calculated so far.  This acts to ensure that the displacement field remains globally smooth even
after many iterations, and is enabled by default.  For many registration operations this is
desirable behaviour: if the images being registered are related by some kind of global
transformation, for example a stretching, scaling or warping as might be encountered when
registering images of two different organs, or growth of a structure over time.  In other cases
this may lead to a poorer registration, particularly if there are localised deformations in the
image pair, examples of such situations might include determining the relative motion of two
structures, or identifying small localised changes over a large image.  The memory term is
controlled by the ``with_memory`` parameter and defaults to ``with_memory = true``.

The cartoon faces we have used for demonstration so far are an excellent example of a situation
where there are localised changes: in the registration of the sad face to the happy face the border
of the face and the eyes are unchanged between the two images and only the mouth changes.  In all
the registration examples so far the global smoothing has caused distortion of these static
structures as well as limiting the accuracy of registration of the mouth area.

The file ``faces_1/sad2happy_no_memory.conf`` disables the memory term and hence the global
smoothing:

.. literalinclude:: /tutorial_files/faces_1/sad2happy_no_memory.conf
   :language: ini

Running the registration with these parameters produces the following result:

+-----------------------------+-----------------------------+
| |sad2happy_no_memory|       | |happy|                     |
+-----------------------------+-----------------------------+
| ``sad2happy_no_memory.png`` | ``happy.png``               |
+-----------------------------+-----------------------------+

.. |sad2happy_no_memory| image:: /tutorial_files/faces_1/sad2happy_no_memory.png

With the global smoothing term disabled the quality of the registration is vastly improved.  The
mouth is cleanly registered with much smaller error, and the border of the face, along with the
eyes remain static through the registration.

*The ``with memory`` parameter provides control over the global smoothing of the displacement
field.  The appropriateness of the global smoothing is problem specific, and for more unusual
problems some experimentation may be required to find the optimal smoothing behaviour.*


Unregisterable Images
=====================

It is important to realise that there are a large class of image pairs for which pFIRE cannot provide a
satisfactory registration.  In order for two images to be registerable they must have all the same
features (topologically speaking they must be homotopic_) such that one can be continuously
deformed into the other.

.. _homotopic: http://mathworld.wolfram.com/Homotopic.html

As a practical example consider the pair of images we registered in the example above.  There, the
sad face can be registered into the happy face by displacing the sides of the mouth upwards and a
smooth displacement field can be used to do this.

In comparison, consider the two images below:

+--------------------+--------------------+
| |grin|             | |sad|              |
+--------------------+--------------------+
| ``grin.png``       | ``sad.png``        |
+--------------------+--------------------+

In this case the grinning face does not have the same features as the sad face.  In order to
distort the sad face into the grinning face the pixels of the mouth would have to be simultaneously
displaced upwards to form the top of the grinning mouth, and downwards to form the bottom of it.
Since the displacement field must be single valued this is not possible, and so the images cannot
be registered.

These images are located in the ``faces_2`` example folder. A sample registration configuration is
provided in ``sad2grin.conf``.  Performing the registration should result in:

+--------------------+--------------------+
| |sad2grin|         | |grin|             |
+--------------------+--------------------+
| ``sad2grin.png``   | ``grin.png``       |
+--------------------+--------------------+

.. |grin| image:: /tutorial_files/faces_2/grin.png

.. |sad2grin| image:: /tutorial_files/faces_2/sad2grin.png

As you can see: *pFIRE will always produce a result. It is entirely up to the user to determine if
two images are suitable to be registered, and to check the results are sane!*
