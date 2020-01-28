.. _em:

Energy Minimization
-------------------

Energy minimization in |Gromacs| can be done using steepest descent,
conjugate gradients, or l-bfgs (limited-memory
Broyden-Fletcher-Goldfarb-Shanno quasi-Newtonian minimizer...we prefer
the abbreviation). EM is just an option of the :ref:`mdrun <gmx mdrun>` program.

Steepest Descent
~~~~~~~~~~~~~~~~

Although steepest descent is certainly not the most efficient algorithm
for searching, it is robust and easy to implement.

We define the vector :math:`\mathbf{r}` as the vector of
all :math:`3N` coordinates. Initially a maximum displacement :math:`h_0`
(*e.g.* 0.01 nm) must be given.

First the forces :math:`\mathbf{F}` and potential energy
are calculated. New positions are calculated by

  .. math:: \mathbf{r}_{n+1} =  \mathbf{r}_n + \frac{\mathbf{F}_n}{\max (|\mathbf{F}_n|)} h_n,
            :label: eqnEMpos

where :math:`h_n` is the maximum displacement and
:math:`\mathbf{F}_n` is the force, or the negative
gradient of the potential :math:`V`. The notation :math:`\max
(|\mathbf{F}_n|)` means the largest scalar force on any
atom. The forces and energy are again computed for the new positions

| If (:math:`V_{n+1} < V_n`) the new positions are accepted and
  :math:`h_{n+1} = 1.2
  h_n`.
| If (:math:`V_{n+1} \geq V_n`) the new positions are rejected and
  :math:`h_n = 0.2 h_n`.

The algorithm stops when either a user-specified number of force
evaluations has been performed (*e.g.* 100), or when the maximum of the
absolute values of the force (gradient) components is smaller than a
specified value :math:`\epsilon`. Since force truncation produces some
noise in the energy evaluation, the stopping criterion should not be
made too tight to avoid endless iterations. A reasonable value for
:math:`\epsilon` can be estimated from the root mean square force
:math:`f` a harmonic oscillator would exhibit at a temperature
:math:`T`. This value is

.. math:: f = 2 \pi \nu \sqrt{ 2mkT},
          :label: eqnEMharmosc

where :math:`\nu` is the oscillator frequency, :math:`m` the (reduced)
mass, and :math:`k` Boltzmann’s constant. For a weak oscillator with a
wave number of 100 cm\ :math:`^{-1}` and a mass of 10 atomic units, at a
temperature of 1 K, :math:`f=7.7` kJ mol\ :math:`^{-1}` nm\ :math:`^{-1}`.
A value for :math:`\epsilon` between 1 and 10 is acceptable.

Conjugate Gradient
~~~~~~~~~~~~~~~~~~

Conjugate gradient is slower than steepest descent in the early stages
of the minimization, but becomes more efficient closer to the energy
minimum. The parameters and stop criterion are the same as for steepest
descent. In |Gromacs| conjugate gradient can not be used with constraints,
including the SETTLE algorithm for water \ :ref:`47 <refMiyamoto92>`, as
this has not been implemented. If water is present it must be of a
flexible model, which can be specified in the :ref:`mdp` file
by ``define = -DFLEXIBLE``.

This is not really a restriction, since the accuracy of conjugate
gradient is only required for minimization prior to a normal-mode
analysis, which cannot be performed with constraints. For most other
purposes steepest descent is efficient enough.

L-BFGS
~~~~~~

The original BFGS algorithm works by successively creating better
approximations of the inverse Hessian matrix, and moving the system to
the currently estimated minimum. The memory requirements for this are
proportional to the square of the number of particles, so it is not
practical for large systems like biomolecules. Instead, we use the
L-BFGS algorithm of Nocedal  \ :ref:`52 <refByrd95a>`,
:ref:`53 <refZhu97a>`, which approximates the inverse Hessian by a fixed number
of corrections from previous steps. This sliding-window technique is
almost as efficient as the original method, but the memory requirements
are much lower - proportional to the number of particles multiplied with
the correction steps. In practice we have found it to converge faster
than conjugate gradients, but due to the correction steps it is not yet
parallelized. It is also noteworthy that switched or shifted
interactions usually improve the convergence, since sharp cut-offs mean
the potential function at the current coordinates is slightly different
from the previous steps used to build the inverse Hessian approximation.
