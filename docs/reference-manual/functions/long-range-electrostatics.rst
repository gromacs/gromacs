.. _lrelstat:   

Long Range Electrostatics
-------------------------

Ewald summation
~~~~~~~~~~~~~~~

The total electrostatic energy of :math:`N` particles and their periodic
images is given by

.. math:: V=\frac{f}{2}\sum_{n_x}\sum_{n_y}
          \sum_{n_{z}*} \sum_{i}^{N} \sum_{j}^{N}
          \frac{q_i q_j}{{\bf r}_{ij,{\bf n}}}.
          :label: eqntotalcoulomb

:math:`(n_x,n_y,n_z)={\bf n}` is the box index vector, and the star
indicates that terms with :math:`i=j` should be omitted when
:math:`(n_x,n_y,n_z)=(0,0,0)`. The distance :math:`{\bf r}_{ij,{\bf n}}`
is the real distance between the charges and not the minimum-image. This
sum is conditionally convergent, but very slow.

Ewald summation was first introduced as a method to calculate long-range
interactions of the periodic images in crystals \ :ref:`105 <refEwald21>`. The idea
is to convert the single slowly-converging sum :eq:`eqn. %s <eqntotalcoulomb>`
into two quickly-converging terms and a constant term:

.. math:: \begin{aligned}
          V &=& V_{\mathrm{dir}} + V_{\mathrm{rec}} + V_{0} \\[0.5ex]
          V_{\mathrm{dir}} &=& \frac{f}{2} \sum_{i,j}^{N}
          \sum_{n_x}\sum_{n_y}
          \sum_{n_{z}*} q_i q_j \frac{\mbox{erfc}(\beta {r}_{ij,{\bf n}} )}{{r}_{ij,{\bf n}}} \\[0.5ex]
          V_{\mathrm{rec}} &=& \frac{f}{2 \pi V} \sum_{i,j}^{N} q_i q_j
          \sum_{m_x}\sum_{m_y}
          \sum_{m_{z}*} \frac{\exp{\left( -(\pi {\bf m}/\beta)^2 + 2 \pi i
                {\bf m} \cdot ({\bf r}_i - {\bf r}_j)\right)}}{{\bf m}^2} \\[0.5ex]
          V_{0} &=& -\frac{f \beta}{\sqrt{\pi}}\sum_{i}^{N} q_i^2,\end{aligned}
          :label: eqntotalcoloumbseparate

where :math:`\beta` is a parameter that determines the relative weight
of the direct and reciprocal sums and :math:`{\bf m}=(m_x,m_y,m_z)`. In
this way we can use a short cut-off (of the order of :math:`1`  nm) in
the direct space sum and a short cut-off in the reciprocal space sum
(*e.g.* 10 wave vectors in each direction). Unfortunately, the
computational cost of the reciprocal part of the sum increases as
:math:`N^2` (or :math:`N^{3/2}` with a slightly better algorithm) and it
is therefore not realistic for use in large systems.

Using Ewald
^^^^^^^^^^^

Don’t use Ewald unless you are absolutely sure this is what you want -
for almost all cases the PME method below will perform much better. If
you still want to employ classical Ewald summation enter this in your
:ref:`mdp` file, if the side of your box is about :math:`3`  nm:

::

    coulombtype     = Ewald
    rvdw            = 0.9
    rlist           = 0.9
    rcoulomb        = 0.9
    fourierspacing  = 0.6
    ewald-rtol      = 1e-5

The ratio of the box dimensions and the fourierspacing parameter
determines the highest magnitude of wave vectors :math:`m_x,m_y,m_z` to
use in each direction. With a 3-nm cubic box this example would use
:math:`11` wave vectors (from :math:`-5` to :math:`5`) in each
direction. The ``ewald-rtol`` parameter is the relative strength of the
electrostatic interaction at the cut-off. Decreasing this gives you a
more accurate direct sum, but a less accurate reciprocal sum.

.. _pme:

PME
~~~

Particle-mesh Ewald is a method proposed by Tom
Darden \ :ref:`14 <refDarden93>` to improve the performance of the reciprocal sum.
Instead of directly summing wave vectors, the charges are assigned to a
grid using interpolation. The implementation in |Gromacs| uses cardinal
B-spline interpolation \ :ref:`15 <refEssmann95>`, which is referred to as
smooth PME (SPME). The grid is then Fourier transformed with a 3D FFT
algorithm and the reciprocal energy term obtained by a single sum over
the grid in k-space.

The potential at the grid points is calculated by inverse
transformation, and by using the interpolation factors we get the forces
on each atom.

The PME algorithm scales as :math:`N \log(N)`, and is substantially
faster than ordinary Ewald summation on medium to large systems. On very
small systems it might still be better to use Ewald to avoid the
overhead in setting up grids and transforms. For the parallelization of
PME see the section on MPMD PME (:ref:`mpmdpme`).

With the Verlet cut-off scheme, the PME direct space potential is
shifted by a constant such that the potential is zero at the cut-off.
This shift is small and since the net system charge is close to zero,
the total shift is very small, unlike in the case of the Lennard-Jones
potential where all shifts add up. We apply the shift anyhow, such that
the potential is the exact integral of the force.

Using PME
^^^^^^^^^

As an example for using Particle-mesh Ewald summation in |Gromacs|,
specify the following lines in your :ref:`mdp` file:

::

    coulombtype     = PME
    rvdw            = 0.9
    rlist           = 0.9
    rcoulomb        = 0.9
    fourierspacing  = 0.12
    pme-order       = 4
    ewald-rtol      = 1e-5

In this case the ``fourierspacing`` parameter determines the
maximum spacing for the FFT grid (i.e. minimum number of grid points),
and ``pme-order`` controls the interpolation order. Using
fourth-order (cubic) interpolation and this spacing should give
electrostatic energies accurate to about :math:`5\cdot10^{-3}`. Since
the Lennard-Jones energies are not this accurate it might even be
possible to increase this spacing slightly.

Pressure scaling works with PME, but be aware of the fact that
anisotropic scaling can introduce artificial ordering in some systems.

P3M-AD
~~~~~~

The Particle-Particle
Particle-Mesh
methods of Hockney & Eastwood can also be applied in |Gromacs| for the
treatment of long range electrostatic
interactions \ :ref:`106 <refHockney81>`. Although the P3M method was the first efficient long-range
electrostatics method for molecular simulation, the smooth PME (SPME)
method has largely replaced P3M as the method of choice in atomistic
simulations. One performance disadvantage of the original P3M method was
that it required 3 3D-FFT back transforms to obtain the forces on the
particles. But this is not required for P3M and the forces can be
derived through analytical differentiation of the potential, as done in
PME. The resulting method is termed P3M-AD. The only remaining
difference between P3M-AD and PME is the optimization of the lattice
Green influence function for error minimization that P3M uses. However,
in 2012 it has been shown that the SPME influence function can be
modified to obtain P3M \ :ref:`107 <refBallenegger2012>`. This means
that the advantage of error minimization in P3M-AD can be used at the
same computational cost and with the same code as PME, just by adding a
few lines to modify the influence function. However, at optimal
parameter setting the effect of error minimization in P3M-AD is less
than 10%. P3M-AD does show large accuracy gains with interlaced (also
known as staggered) grids, but that is not supported in |Gromacs| (yet).

P3M is used in |Gromacs| with exactly the same options as used with PME by
selecting the electrostatics type:

::

    coulombtype     = P3M-AD

Optimizing Fourier transforms and PME calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to optimize the parameters for calculation of
electrostatic interaction such as PME grid dimensions and cut-off radii.
This is particularly relevant to do before launching long production
runs.

:ref:`gmx mdrun` will automatically do a lot of PME
optimization, and |Gromacs| also includes a special tool,
:ref:`gmx tune_pme`, which automates the process of selecting
the optimal number of PME-only ranks.
