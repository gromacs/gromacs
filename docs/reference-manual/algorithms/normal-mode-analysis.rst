Normal-Mode Analysis
--------------------

Normal-mode analysis \ :ref:`54 <refLevitt83>`\ :ref:`56 <refBBrooks83b>`
can be performed using |Gromacs|, by diagonalization of the
mass-weighted Hessian :math:`H`:

.. math:: \begin{aligned}
          R^T M^{-1/2} H M^{-1/2} R   &=& \mbox{diag}(\lambda_1,\ldots,\lambda_{3N})
          \\
          \lambda_i &=& (2 \pi \omega_i)^2\end{aligned}
          :label: eqnNMA

where :math:`M` contains the atomic masses, :math:`R` is a matrix that
contains the eigenvectors as columns, :math:`\lambda_i` are the
eigenvalues and :math:`\omega_i` are the corresponding frequencies.

First the Hessian matrix, which is a :math:`3N \times 3N` matrix where
:math:`N` is the number of atoms, needs to be calculated:

.. math:: \begin{aligned}
          H_{ij}  &=&     \frac{\partial^2 V}{\partial x_i \partial x_j}\end{aligned}
          :label: eqnNMAhessian

where :math:`x_i` and :math:`x_j` denote the atomic x, y or z
coordinates. In practice, this equation is not used, but the Hessian is
calculated numerically from the force as:

.. math:: \begin{aligned}
          H_{ij} &=& -
            \frac{f_i({\bf x}+h{\bf e}_j) - f_i({\bf x}-h{\bf e}_j)}{2h}
          \\
          f_i     &=& - \frac{\partial V}{\partial x_i}\end{aligned}
          :label: eqnNMAhessianfromforce

where :math:`{\bf e}_j` is the unit vector in direction :math:`j`. It
should be noted that for a usual normal-mode calculation, it is
necessary to completely minimize the energy prior to computation of the
Hessian. The tolerance required depends on the type of system, but a
rough indication is 0.001 kJ mol\ :math:`^{-1}`. Minimization should be
done with conjugate gradients or L-BFGS in double precision.

A number of |Gromacs| programs are involved in these calculations. First,
the energy should be minimized using :ref:`mdrun <gmx mdrun>`. Then,
:ref:`mdrun <gmx mdrun>` computes the Hessian. **Note** that for generating
the run input file, one should use the minimized conformation from the
full precision trajectory file, as the structure file is not accurate
enough. :ref:`gmx nmeig` does the
diagonalization and the sorting of the normal modes according to their
frequencies. Both :ref:`mdrun <gmx mdrun>` and :ref:`gmx nmeig` should be run in double precision.
The normal modes can be analyzed with the program :ref:`gmx anaeig`. Ensembles
of structures at any temperature and for any subset of normal modes can
be generated with :ref:`gmx nmens`. An overview of normal-mode analysis and the
related principal component analysis (see sec. :ref:`covanal`) can be
found in \ :ref:`57 <refHayward95b>`.
