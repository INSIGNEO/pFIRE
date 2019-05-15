===================
The pFIRE Algorithm
===================

The basic equation that pFIRE attempts to solve relates the fixed image :math:`\vec{f}` to the
moved image :math:`\vec{m}` via the displacement map :math:`\vec{a}`.  This can be written in one
of two forms, depending on if we are "pushing" the moved image onto the fixed image or pulling the
fixed onto the moved image.

.. math:: 

  \vec{f}(\vec{x}) = \vec{m}(\vec{x}+\vec{a}(\vec{x}))

  \vec{m}(\vec{x}) = \vec{f}(\vec{x}-\vec{a}(\vec{x}))

Taylor expanding these expressions and truncating to first order allows us to separate the
displacement vector from the source image by introducing a gradient term:

.. math::

  \vec{f}(\vec{x}) \simeq \vec{m}(\vec{x}) + \nabla\vec{m}(\vec{x}) \cdot \vec{a}(\vec{x})

  \vec{m}(\vec{x}) \simeq \vec{f}(\vec{x}) - \nabla\vec{f}(\vec{x}) \cdot \vec{a}(\vec{x})

Note that these expressions are valid only in the case that the displacement
:math:`\vec{a}(\vec{x})` is small, however, the scheme will be applied iteratively to converge upon
the correct result.

Combining these two expressions yields a single registration equation to be solved for
:math:`\vec{a}(\vec{x})` which includes the information from both images:

.. math::
 
  \vec{f} - \vec{m} = \frac{1}{2}\nabla[\vec{f}(\vec{x}) + \vec{m}(\vec{x})]\cdot\vec{a}(\vec{x}).

The above registration equation has the images and map in the form on continuous functions,
however, when solving numerically we have images composed of pixels and a map represented as a
series of node points.  Information between these points can then be found by (linear)
interpolation.  Typically there will also be far fewer map nodes than image pixels both for
practical reasons of computation time and mathematical reasons of resistance to noise and simply
having sufficient information to fully constrain the problem.

The net result of this is that we require some kind of interpolation method to map contributions
from image pixels to their corresponding map nodes.  pFIRE uses linear interpolation, with each
pixel contributing to its nearest 4 (or 8 in 3D) map nodes, with linear interpolation coefficients
given by

.. math::

  \phi_{\mu\nu\omega}^{ijk} = \begin{cases}
    \prod_{n=1,2,3} (1 - | x^n_{ijk} - x^n_{\mu\nu\omega}| ) &\mathrm{if}\ 
      \forall\ (x^n_{ijk} - x^n_{\mu\nu\omega}) < 1\\
    0 &\mathrm{otherwise}
  \end{cases}

Here upper latin indices represent the image pixel indices and the greek indices represent the map
node indices. The object :math:`x^n_{ijk}` is the position in spatial dimension :math:`n` of the
pixel at position (i,j,k) in the image, and similarly :math:`x^n_{\mu\nu\omega}` describes the
position of a map node.

This can then be used, along with a gradient information tensor :math:`I` and a displacement (map)
tensor :math:`A`, to construct a discrete form of the registration equation (using Einstein
summation convention):

.. math::

  (f - m)_{xyz} = I^{xyz}_{\alpha ijk} \phi^{ijk}_{\mu\nu\omega} A^{\alpha \mu\nu\omega}

where the index :math:`\alpha` refers to the four degrees of freedom (:math:`x, y, z, s`), the
latin index triplets (:math:`x y z` and :math:`i j k`) refer to image pixel indices and the greek
index triplet (:math:`\mu \nu \omega`) refers to map node indices.

The form of the gradient information tensor :math:`I` is diagonal, with elements given by

.. math::

  I^{xyz}_{\alpha ijk} = \begin{cases}
   \frac{1}{2}(f + m)_{xyz} & \text{if}\ x=i, y=j, z=k, \alpha=0\\
   \frac{1}{2}\nabla_\alpha (f + m)_{xyz} & \text{if}\ x=i, y=j, z=k, \alpha={1,2,3}\\
   0 & \text{otherwise}
  \end{cases}

where :math:`\nabla_\alpha` indicates differentiation with respect to :math:`x, y\ \text{or}\ z`
for :math:`\alpha=1,2,3` respectively. Note that this is essentially a concatenation of four
diagonal tensors, and in the subsequent steps will become a concatentation of four diagonal
matrices.


The complexity of implementation may be eased by electing to explicitly flatten the indices
:math:`x,y,z` and :math:`\mu,\nu,\omega`, using either row-major or column-major indexing,
depending on the implementation requirements.

This allows the problem to be expressed in the form of a matrix equation:

.. math::
  
  \begin{align}
    (\bar{F} - \bar{M}) &= \bar{\bar{T}}\bar{A} = \bar{\bar{I}}\bar{\bar{\Phi}}\bar{A} \\
    &= \begin{bmatrix} \bar{\bar{I}}_s \\ \bar{\bar{I}}_x \\ \bar{\bar{I}}_y \\ \bar{\bar{I}}_z \end{bmatrix}
    \begin{bmatrix} \bar{\bar{\phi}} \\ \bar{\bar{\phi}} \\ \bar{\bar{\phi}} \\ \bar{\bar{\phi}} \end{bmatrix}
    \begin{bmatrix} \bar{A}_s \\ \bar{A}_x \\ \bar{A}_y \\ \bar{A}_z \end{bmatrix}
  \end{align}

where :math:`\bar{F}` and :math:`\bar{M}` are column vectors containing the fixed and moved image
intensities, :math:`\bar{\bar{I}}_\alpha` are the square, diagonal submatrices for each dimension
of the gradient information, :math:`\bar{\bar{\phi}}` is the same submatrix for each dimension, and
:math:`A_\alpha` are column matrices containing the nodal displacements for each dimension of the
map.

The matrix :math:`T` is defined here as the product of the matrices :math:`\bar{\bar{I}}` and
:math:`\bar{\bar{\Phi}}`, and it is, in general, non-square.  We therefore create a soluable matrix
equation from this by multiplying both sides by the transpose :math:`\bar{\bar{T}}^t`,

.. math::

  \bar{\bar{T}}^t(\bar{F} - \bar{M}) = \bar{\bar{T}}^t\bar{\bar{T}}\bar{A}

yielding an equation with square matrix :math:`\bar{\bar{T}}^t\bar{\bar{T}}` which may be solved by
standard numerical methods.

Unfortunately :math:`\bar{\bar{T}}^t\bar{\bar{T}}` is typically singular and so we must apply some
kind of regularisation to make the problem soluable.  In this case we use the method of Tikhonov
regularisation, which includes an additional smoothing term in the problem. In this case the
Laplacian matrix is chosen, adding minimization of the second derivative of the map as the extra
constraint (:math:`\bar{\bar{L}}\bar{A} = 0` and the problem becomes,

.. math::

 \begin{bmatrix}\bar{\bar{T}}^t & \lambda^\frac{1}{2}\bar{\bar{L}}^t \end{bmatrix}
 \begin{bmatrix}(\bar{F} - \bar{M}) \\ \bar{\bar{0}}\end{bmatrix}
  = \begin{bmatrix}\bar{\bar{T}}^t\bar{\bar{T}} + \lambda\bar{\bar{L}}^t\bar{\bar{L}}\end{bmatrix}\bar{A}.

This simplifies to

.. math::

 \bar{\bar{T}}^t(\bar{F} - \bar{M})
  = \begin{bmatrix}\bar{\bar{T}}^t\bar{\bar{T}} + \lambda\bar{\bar{L}}^t\bar{\bar{L}}\end{bmatrix}\bar{A}

where :math:`\lambda` is a free parameter which is chosen to minimize the condition number of the
regularised matrix. Adding this additional Laplacian constraint ensures that the smoothest solution
to the problem is chosen.

Finally, because the algorithm is iterative, we calculate a new mapping at each step and add it to
the previously accumulated displacements to refine the overall solution.  In many cases, we may want to
constrain the overall displacement to be maximally smooth, thus our additional constraint becomes

.. math::

  \bar{\bar{L}}(\bar{A} + \bar{A}_p) = 0

or

.. math::

  \bar{\bar{L}}(\bar{A} = - \bar{A}_p)


So the equation becomes

.. math::

 \begin{bmatrix}\bar{\bar{T}}^t & \lambda^\frac{1}{2}\bar{\bar{L}}^t \end{bmatrix}
 \begin{bmatrix}(\bar{F} - \bar{M}) \\ -\lambda^\frac{1}{2}\bar{\bar{L}}^t\bar{A}_p\end{bmatrix}
  = \begin{bmatrix}\bar{\bar{T}}^t\bar{\bar{T}} + \lambda\bar{\bar{L}}^t\bar{\bar{L}}\end{bmatrix}\bar{A}.

This simplifies to

.. math::

 \bar{\bar{T}}^t(\bar{F} - \bar{M}) - \lambda\bar{\bar{L}}^t\bar{\bar{L}}\bar{A}_p
  = \begin{bmatrix}\bar{\bar{T}}^t\bar{\bar{T}} + \lambda\bar{\bar{L}}^t\bar{\bar{L}}\end{bmatrix}\bar{A}

This form of the equation "remembers" the value of the map from the previous iteration and attempts
to enforce global smoothing on the final result and is referred to as the "memory term".  This
method of smoothing is useful for registering images where the displacement is expected to be
continuous and smooth, for example in the case of registration of multimodal images of the same
structure.  In constrast, this option should be disabled in the case that image features are
expected to move relative to each other, for example in cell-tracking applications.
