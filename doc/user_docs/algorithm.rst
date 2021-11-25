===================
The pFIRE Algorithm
===================

---------------------
Algorithm Description
---------------------

pFIRE operates by iteratively solving a matrix equation to find updates to the deformation map that
minimizes the difference between the fixed and the moved image. After each iteration the original
moved image is warped using the updated map to produce a new moved image that is closer to the
fixed image.  The original moved image is used each time to minimize the accumulation of errors
from warping the same image multiple times.

Initially the displacement map that is used is coarse, having only 3 nodes per dimension, but as
the registration converges at each resolution the map is then refined, doubling the number of nodes
per dimension until the desired resolution is reached.

The Registration Equation
-------------------------

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
:math:`x,y,z` and :math:`\mu,\nu,\omega`. In the case of pFIRE indexing is performed such that the
x index varies fastest - this is the same indexing scheme as used by PETSc:

.. math::

  n(x, y, z) = (z * Y + y)*X + x

where :math:`(X,Y,Z)` is the size of the mesh.

This allows the problem to be expressed in the form of a matrix equation:

.. math::

  \begin{align}
    (\bar{F} - \bar{M}) &= \mathbf{{T}}\bar{A} = \mathbf{{G}}\mathbf{{\Phi}}\bar{A} \\
    &= \begin{bmatrix} \mathbf{{G}}_s \\ \mathbf{{G}}_x \\ \mathbf{{G}}_y \\ \mathbf{{G}}_z \end{bmatrix}
    \begin{bmatrix} \mathbf{{\phi}} \\ \mathbf{{\phi}} \\ \mathbf{{\phi}} \\ \mathbf{{\phi}} \end{bmatrix}
    \begin{bmatrix} \bar{A}_s \\ \bar{A}_x \\ \bar{A}_y \\ \bar{A}_z \end{bmatrix}
  \end{align}

where :math:`\bar{F}` and :math:`\bar{M}` are column vectors containing the fixed and moved image
intensities, :math:`\mathbf{{G}}_\alpha` are the square, diagonal submatrices for each dimension
of the gradient information, :math:`\mathbf{{\phi}}` is the same submatrix for each dimension, and
:math:`A_\alpha` are column matrices containing the nodal displacements for each dimension of the
map.

The matrix :math:`T` is defined here as the product of the matrices :math:`\mathbf{{G}}` and
:math:`\mathbf{{\Phi}}`, and it is, in general, non-square.  We therefore create a solvable matrix
equation from this by multiplying both sides by the transpose :math:`\mathbf{{T}}^t`,

.. math::

  \mathbf{{T}}^t(\bar{F} - \bar{M}) = \mathbf{{T}}^t\mathbf{{T}}\bar{A}

yielding an equation with square matrix :math:`\mathbf{{T}}^t\mathbf{{T}}` which may be solved by
standard numerical methods.

Unfortunately :math:`\mathbf{{T}}^t\mathbf{{T}}` is typically singular and so we must apply some
kind of regularisation to make the problem solvable.  In this case we use the method of Tikhonov
regularisation, which includes an additional smoothing term in the problem. In this case the
Laplacian matrix is chosen, adding minimization of the second derivative of the map as the extra
constraint (:math:`\mathbf{{L}}\bar{A} = 0` and the problem becomes,

.. math::

 \begin{bmatrix}\mathbf{{T}}^t & \lambda^\frac{1}{2}\mathbf{{L}}^t \end{bmatrix}
 \begin{bmatrix}(\bar{F} - \bar{M}) \\ \mathbf{{0}}\end{bmatrix}
  = \begin{bmatrix}\mathbf{{T}}^t\mathbf{{T}} + \lambda\mathbf{{L}}^t\mathbf{{L}}\end{bmatrix}\bar{A}.

This simplifies to

.. math::

 \mathbf{{T}}^t(\bar{F} - \bar{M})
  = \begin{bmatrix}\mathbf{{T}}^t\mathbf{{T}} + \lambda\mathbf{{L}}^t\mathbf{{L}}\end{bmatrix}\bar{A}

where :math:`\lambda` is a free parameter which is chosen to minimize the condition number of the
regularised matrix. Adding this additional Laplacian constraint ensures that the smoothest solution
to the problem is chosen.

Finally, because the algorithm is iterative, we calculate a new mapping at each step and add it to
the previously accumulated displacements to refine the overall solution.  In many cases, we may want to
constrain the overall displacement to be maximally smooth, thus our additional constraint becomes

.. math::

  \mathbf{{L}}(\bar{A} + \bar{A}_p) = 0

or

.. math::

  \mathbf{{L}}(\bar{A} = - \bar{A}_p)


So the equation becomes

.. math::

 \begin{bmatrix}\mathbf{{T}}^t & \lambda^\frac{1}{2}\mathbf{{L}}^t \end{bmatrix}
 \begin{bmatrix}(\bar{F} - \bar{M}) \\ -\lambda^\frac{1}{2}\mathbf{{L}}^t\bar{A}_p\end{bmatrix}
  = \begin{bmatrix}\mathbf{{T}}^t\mathbf{{T}} + \lambda\mathbf{{L}}^t\mathbf{{L}}\end{bmatrix}\bar{A}.

This simplifies to

.. math::

 \mathbf{{T}}^t(\bar{F} - \bar{M}) - \lambda\mathbf{{L}}^t\mathbf{{L}}\bar{A}_p
  = \begin{bmatrix}\mathbf{{T}}^t\mathbf{{T}} + \lambda\mathbf{{L}}^t\mathbf{{L}}\end{bmatrix}\bar{A}

This form of the equation "remembers" the value of the map from the previous iteration and attempts
to enforce global smoothing on the final result and is referred to as the "memory term".  This
method of smoothing is useful for registering images where the displacement is expected to be
continuous and smooth, for example in the case of registration of multimodal images of the same
structure.  In constrast, this option should be disabled in the case that image features are
expected to move relative to each other, for example in cell-tracking applications.

----------------------
Implementation Details
----------------------

Calculating :math:`\mathbf{{T}}^t\mathbf{{T}}`
----------------------------------------------

The naive implementation of registration equation requires explicitly constructing the matrix
:math:`\mathbf{{T}}` before calculating the product :math:`\mathbf{{T}}^t\mathbf{{T}}`.  This is
not an efficient use of memory, however, especially since the matrix :math:`\mathbf{{T}}` is much
larger than the final product and is only an intermediate value in computing
:math:`\mathbf{{T}}^t(\bar{F} - \bar{M})` and :math:`\begin{bmatrix}\mathbf{{T}}^t\mathbf{{T}} +
\lambda\mathbf{{L}}^t\mathbf{{L}}\end{bmatrix}`.

Because of this, pFIRE computes the matrix :math:`\mathbf{{T}}^t\mathbf{{T}}` and the vector
:math:`\mathbf{{T}}^t(\bar{F} - \bar{M})` directly.  The structure of the final matrix and vector
are distributed between the various ranks such that each rank computes an equal number of
components of the final data structure.  For each matrix or vector element, the location of all the
required image or gradient pixels is determiend and then required values that are not local to the
rank are communicated using MPI_Alltoall with all ranks communicating at once.

Implementation of the computation can be made more efficient by understanding the structure of the
:math:`\mathbf{{G}}` and :math:`\mathbf{{\Phi}}` matrices, as prior knowledge of the zero-patterns
in these matrices can make calculation of the final matrix :math:`\mathbf{{T}}^t\mathbf{{T}}` much
more efficient.  Further, since the :math:`\mathbf{{G}}` matrix is diagonal we know that the zero
pattern of its product with any matrix will match that of the other matrix, therefore we need only
consider the zero-pattern of the interpolation matrix :math:`\mathbf{{\Phi}}`.

Revisiting the definition of :math:`\mathbf{\phi}` (for flattened indices) and a single dimension:

.. math::

  \phi_{\mu}^{i} = \begin{cases}
    \prod_{n=1,2,3} (1 - | x^n_{i} - x^n_{\mu}| ) &\mathrm{if}\
      \forall\ (x^n_{i} - x^n_{\mu}) < 1\\
    0 &\mathrm{otherwise}
  \end{cases}

with :math:`\mu` the map indices and :math:`i` the image indices. So when calculating

.. math::

  \mathbf{{T}}^t\mathbf{{T}} = \mathbf{\phi}^\mu_i\mathbf{G}^i_j\mathbf{G}^j_k\mathbf{\phi}^k_\nu

the :math:`\mathbf{G}` matrices function as Kronecker deltas for the purposes of discerning the
zero-pattern, and so we can determine a zero-pattern matrix :math:`\mathbf{Z}` for the matrix product
:math:`\mathbf{\phi}^t\mathbf{\phi}`

.. math::

  \mathbf{Z} = \sum_i
  \begin{cases}
      1 &\mathrm{if}\ \forall\ (x^n_{i} - x^n_{\mu}) < 1\\
      0 &\mathrm{otherwise}
    \end{cases}
    \times
    \begin{cases}
     1 &\mathrm{if}\ \forall\ (x^n_{i} - x^n_{\nu}) < 1\\
     0 &\mathrm{otherwise}
    \end{cases}.

From this we see that there will only be a non-zero entry in the map if the image pixel at point
:math:`i` is a contributor to the map node at both points :math:`\mu` and :math:`\nu`.  Because we
are using nearest neighbor interpolation, this means there are two classes of non-zero entry in
:math:`\mathbf{\phi}`: diagonal entries, for which :math:`\mu = \nu`; and off diagonal entries, for
which :math:`\mu = \nu \pm \{1,U,U*V\}` (where :math:`U` and :math:`U*V` are the two multiplicative
factors used in flattening the indexing). This gives a total of 9 non-zero entries per row.  If the
domain splitting is carefully chosen, this can lead to an optimal communication pattern where
information need only be exchanged with nearest neighbour ranks, rather than requiring all-to-all
communication. In practice this optimal scheme is not implemented by pFIRE since it requires fixing
the problem size per processor. Instead dynamic communication patterns are used to communicate only
the required data regardless of its location.

Improving Matrix Conditioning
-----------------------------

The matrix as constructed, even with the addition of the laplacian is suboptimally conditioned.  In
order to improve the conditioning further, a preconditioning step is applied to the calculated
matrix :math:`\mathbf{{T}}^t\mathbf{{T}}`, predicated on the fact that the matrix is symmetric and
positive definite, and that the intensity values calculated for the map will not be used to update
the map, but are relevant purely for the correct calculation of the displacement components.

The average of the diagonal elements is calculated for the displacement components and for the
intensity components, and the intensity components rescaled by multiplying with a diagonal matrix
with appropriate values on the diagonal, such that the average of the diagonal intensity elements
is equal to the average for the displacement elements.

Calculating Lambda
------------------

Prior to calculating :math:`\lambda` a "premultiplication factor" :math:`\chi` is calculated and
applied to rescale the laplacian matrix to make the average of the diagonal elements of
:math:`\mathbf{{L}}^t\mathbf{{L}}` equal to the average of the diagonal elements of
:math:`\mathbf{{T}}^t\mathbf{{T}}`.

.. math::

  \chi = \Sum_i [\mathbf{{T}}^t\mathbf{{T}}]_{i,i} / \Sum_i [\mathbf{{L}}^t\mathbf{{L}}]_{i,i}

The optimal value of :math:`\lambda` is then found by minimising the condition number of the matrix

.. math::                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                               
  \mathbf{{T}}^t\mathbf{{T}} + \chi\lambda \mathbf{{L}}^t\mathbf{{L}}].

In general this could be done by various optimisation strategies, however, in pFIRE the approach
taken is to minimize the required number of condition number calculations as they are
computationall expensive. It has been empirically observed that the behaviour of the condition
number with increasing :math:`\lambda` is to initially decrease to a wide minimum region before
increasing again.  Therefore, a reasonable value for lambda can be found by evaluating the
condition number at three points, one expected to be to well below the optimum value, one expected
to be close to it and one expected to be well above. Provided the middle value produces a lower
condition number than the edges, a quadratic fit to these points will then provide a reasonable
estimation of the optimum lambda.  In the case that the middle point is not lower, the gradient of
the fitted line will inform whether to increase or decrease the values of the search points and try
again.

Solving the Registration Equation
---------------------------------

The registration equation is solved using the Krylov subspace solver routines provided by PETsC.

Warping the Image
-----------------

The image is warped by the "pulling method".  For a given pixel in the target image, the map value
is interpolated to the pixel location and used to determine the location of the source pixel in the
source image.  The value at this point is then calculated, using MPI_AllToAll to communicate the
data from other ranks as needed, and inserted into the target image.
