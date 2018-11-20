.. _covanal:

Covariance analysis
-------------------

Covariance analysis, also called principal component analysis or
essential dynamics :ref:`169 <refAmadei93>`\ , can find
correlated motions. It uses the covariance matrix :math:`C` of the
atomic coordinates:

.. math:: C_{ij} = \left \langle 
          M_{ii}^{\frac{1}{2}} (x_i - \langle x_i \rangle)
          M_{jj}^{\frac{1}{2}}  (x_j - \langle x_j \rangle)
          \right \rangle
          :label: eqncovmatrixcoord

where :math:`M` is a diagonal matrix containing the masses of the atoms
(mass-weighted analysis) or the unit matrix (non-mass weighted
analysis). :math:`C` is a symmetric :math:`3N \times 3N` matrix, which
can be diagonalized with an orthonormal transformation matrix :math:`R`:

.. math:: R^T C R = \mbox{diag}(\lambda_1,\lambda_2,\ldots,\lambda_{3N})
          ~~~~\mbox{where}~~\lambda_1 \geq \lambda_2 \geq \ldots \geq \lambda_{3N}
          :label: eqnorthnormtransformmatrix

The columns of :math:`R` are the eigenvectors, also called principal or
essential modes. :math:`R` defines a transformation to a new coordinate
system. The trajectory can be projected on the principal modes to give
the principal components :math:`p_i(t)`:

.. math:: {\bf p}(t) = R^T M^{\frac{1}{2}} ({\bf x}(t) - \langle {\bf x} \rangle)
          :label: eqnprinccomponents

The eigenvalue :math:`\lambda_i` is the mean square fluctuation of
principal component :math:`i`. The first few principal modes often
describe collective, global motions in the system. The trajectory can be
filtered along one (or more) principal modes. For one principal mode
:math:`i` this goes as follows:

.. math:: {\bf x}^f(t) =
          \langle {\bf x} \rangle + M^{-\frac{1}{2}} R_{ * i} \, p_i(t)
          :label: eqnprincmodei

When the analysis is performed on a macromolecule, one often wants to
remove the overall rotation and translation to look at the internal
motion only. This can be achieved by least square fitting to a reference
structure. Care has to be taken that the reference structure is
representative for the ensemble, since the choice of reference structure
influences the covariance matrix.

One should always check if the principal modes are well defined. If the
first principal component resembles a half cosine and the second
resembles a full cosine, you might be filtering noise (see below). A
good way to check the relevance of the first few principal modes is to
calculate the overlap of the sampling between the first and second half
of the simulation. **Note** that this can only be done when the same
reference structure is used for the two halves.

A good measure for the overlap has been defined in \ :ref:`170 <refHess2002b>`. The
elements of the covariance matrix are proportional to the square of the
displacement, so we need to take the square root of the matrix to
examine the extent of sampling. The square root can be calculated from
the eigenvalues :math:`\lambda_i` and the eigenvectors, which are the
columns of the rotation matrix :math:`R`. For a symmetric and
diagonally-dominant matrix :math:`A` of size :math:`3N \times 3N` the
square root can be calculated as:

.. math:: A^\frac{1}{2} = 
          R \, \mbox{diag}(\lambda_1^\frac{1}{2},\lambda_2^\frac{1}{2},\ldots,\lambda_{3N}^\frac{1}{2}) \, R^T
          :label: eqnmatrixsquareroot

It can be verified easily that the product of this matrix with itself
gives :math:`A`. Now we can define a difference :math:`d` between
covariance matrices :math:`A` and :math:`B` as follows:

.. math:: \begin{aligned}
          d(A,B) & = & \sqrt{\mbox{tr}\left(\left(A^\frac{1}{2} - B^\frac{1}{2}\right)^2\right)
          }
          \\ & = &
          \sqrt{\mbox{tr}\left(A + B - 2 A^\frac{1}{2} B^\frac{1}{2}\right)}
          \\ & = &
          \left( \sum_{i=1}^N \left( \lambda_i^A + \lambda_i^B \right)
          - 2 \sum_{i=1}^N \sum_{j=1}^N \sqrt{\lambda_i^A \lambda_j^B}
          \left(R_i^A \cdot R_j^B\right)^2 \right)^\frac{1}{2}\end{aligned}
          :label: eqnmatrixdiff

where tr is the trace of a matrix. We can now define the overlap
:math:`s` as:

.. math:: s(A,B) = 1 - \frac{d(A,B)}{\sqrt{\mbox{tr}A + \mbox{tr} B}}
          :label: eqnmatrixoverlap

The overlap is 1 if and only if matrices :math:`A` and :math:`B` are
identical. It is 0 when the sampled subspaces are completely orthogonal.

A commonly-used measure is the subspace overlap of the first few
eigenvectors of covariance matrices. The overlap of the subspace spanned
by :math:`m` orthonormal vectors :math:`{\bf w}_1,\ldots,{\bf w}_m` with
a reference subspace spanned by :math:`n` orthonormal vectors
:math:`{\bf v}_1,\ldots,{\bf v}_n` can be quantified as follows:

.. math:: \mbox{overlap}({\bf v},{\bf w}) =
          \frac{1}{n} \sum_{i=1}^n \sum_{j=1}^m ({\bf v}_i \cdot {\bf w}_j)^2
          :label: eqnsubspaceoverlap

The overlap will increase with increasing :math:`m` and will be 1 when
set :math:`{\bf v}` is a subspace of set :math:`{\bf w}`. The
disadvantage of this method is that it does not take the eigenvalues
into account. All eigenvectors are weighted equally, and when degenerate
subspaces are present (equal eigenvalues), the calculated overlap will
be too low.

Another useful check is the cosine content. It has been proven that the
the principal components of random diffusion are cosines with the number
of periods equal to half the principal component
index \ :ref:`170 <refHess2002b>`, :ref:`171 <refHess2000>`.
The eigenvalues are proportional to the index to the power
:math:`-2`. The cosine content is defined as:

.. math:: \frac{2}{T}
          \left( \int_0^T \cos\left(\frac{i \pi t}{T}\right) \, p_i(t) \mbox{d} t \right)^2
          \left( \int_0^T p_i^2(t) \mbox{d} t \right)^{-1}
          :label: eqneigenvaluecosine

When the cosine content of the first few principal components is close
to 1, the largest fluctuations are not connected with the potential, but
with random diffusion.

The covariance matrix is built and diagonalized by
:ref:`gmx covar <gmx covar>`. The principal components and
overlap (and many more things) can be plotted and analyzed with
:ref:`gmx anaeig <gmx anaeig>`. The cosine
content can be calculated with
:ref:`gmx analyze <gmx analyze>`.
