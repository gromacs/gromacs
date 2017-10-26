Definitions and Units
=====================

Notation
--------

The following conventions for mathematical typesetting are used
throughout this document:

We define the *lowercase* subscripts :math:`i`, :math:`j`, :math:`k` and
:math:`l` to denote particles: :math:`{\mbox{\boldmath ${r}$}}_i` is the
*position vector* of particle :math:`i`, and using this notation:

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}_ij	=	{\mbox{\boldmath ${r}$}}_j-{{\mbox{\boldmath ${r}$}}_i}\\
   {r_{ij}}=	| {\mbox{\boldmath ${r}$}}_ij | \end{aligned}

The force on particle :math:`i` is denoted by
:math:`{\mbox{\boldmath ${F}$}}_i` and

.. math:: {\mbox{\boldmath ${F}$}}_{ij} = \mbox{force on $i$ exerted by $j$}

Please note that we changed notation as of version 2.0 to
:math:`{\mbox{\boldmath ${r}$}}_ij={\mbox{\boldmath ${r}$}}_j-{{\mbox{\boldmath ${r}$}}_i}`
since this is the notation commonly used. If you encounter an error, let
us know.

MD units
--------

GROMACS uses a consistent set of units that produce values in the
vicinity of unity for most relevant molecular quantities. Let us call
them *MD units*. The basic units in this system are nm, ps, K, electron
charge (e) and atomic mass unit (u), see Table
The values used in GROMACS are
taken from the CODATA Internationally recommended 2010 values of
fundamental physical constants (see `NIST homepage <http://nist.gov>`__). Consistent
with these units are a set of derived units, given in
Table
The **electric conversion factor**
:math:`f=\frac{1}{4 \pi \varepsilon_o}={138.935\,458}`
:math:`\mathrm{kJ}~\mathrm{mol}^{-1}\mathrm{nm}~\mathrm{ e}^{-2}`.
It relates the mechanical quantities to the electrical quantities as in

.. math:: V = f \frac{q^2}{r} \mbox{\ \ or\ \ } F = f \frac{q^2}{r^2}

Electric potentials :math:`\Phi` and electric fields
:math:`{\mbox{\boldmath ${E}$}}` are intermediate quantities in the
calculation of energies and forces. They do not occur inside GROMACS. If
they are used in evaluations, there is a choice of equations and related
units. We strongly recommend following the usual practice of including
the factor :math:`f` in expressions that evaluate :math:`\Phi` and
:math:`{\mbox{\boldmath ${E}$}}`:

.. math::

   \begin{aligned}
   \Phi({\mbox{\boldmath ${r}$}}) = f \sum_j \frac{q_j}{| {\mbox{\boldmath ${r}$}}-{\mbox{\boldmath ${r}$}}_j | } 	\\
   {\mbox{\boldmath ${E}$}}({\mbox{\boldmath ${r}$}}) = f \sum_j q_j \frac{({\mbox{\boldmath ${r}$}}-{\mbox{\boldmath ${r}$}}_j)}{| {\mbox{\boldmath ${r}$}}-{\mbox{\boldmath ${r}$}}_j| ^3}\end{aligned}

With these definitions, :math:`q\Phi` is an energy and
:math:`q{\mbox{\boldmath ${E}$}}` is a force. The units are those given
in Table
about 10 mV for potential.
Thus, the potential of an electronic charge at a distance of 1 nm equals
:math:`f \approx 140` units :math:`\approx 1.4` V.
(exact value: :math:`1.439\,964\,5` V)

**Note** that these units are mutually consistent; changing any of the
units is likely to produce inconsistencies and is therefore *strongly
discouraged*! In particular: if Å are used instead of nm, the unit of
time changes to 0.1 ps. If :math:`\mathrm{kcal}~\mathrm{mol}^{-1}` (= 4.184
:math:`\mathrm{kJ~mol}^{-1}`) is used instead of :math:`\mathrm{kJ~mol}^{-1}` for energy,
the unit of time becomes 0.488882 ps and the unit of temperature changes
to 4.184 K. But in both cases all electrical energies go wrong, because
they will still be computed in :math:`\mathrm{kJ~mol}^{-1}`, expecting nm as
the unit of length. Although careful rescaling of charges may still
yield consistency, it is clear that such confusions must be rigidly
avoided.

In terms of the MD units, the usual physical constants take on different
values (see Table
). All quantities are per
mol rather than per molecule. There is no distinction between
Boltzmann’s constant :math:`k` and the gas constant :math:`R`: their
value is :math:`0.008\,314\,462\,1\mathrm{kJ~mol}^{-1} \mathrm{K}^{-1}`.

Reduced units
-------------

When simulating Lennard-Jones (LJ) systems, it might be advantageous to
use reduced units (*i.e.*, setting
:math:`\epsilon_{ii}=\sigma_{ii}=m_i=k_B=1` for one type of atoms). This
is possible. When specifying the input in reduced units, the output will
also be in reduced units. The one exception is the *temperature*, which
is expressed in :math:`0.008\,314\,462\,1` reduced units. This is a
consequence of using Boltzmann’s constant in the evaluation of
temperature in the code. Thus not :math:`T`, but :math:`k_BT`, is the
reduced temperature. A GROMACS temperature :math:`T=1` means a reduced
temperature of :math:`0.008\ldots` units; if a reduced temperature of 1
is required, the GROMACS temperature should be :math:`120.272\,36`.

In Table
quantities are given for LJ
potentials:

.. math:: V_{LJ} = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]

Mixed or Double precision
-------------------------

|Gromacs| can be compiled in either mixed
or double
precision.
Documentation of previous |Gromacs| versions referred to *single
precision*, but the implementation has made selective use of double
precision for many years. Using single precision for all variables would
lead to a significant reduction in accuracy. Although in *mixed
precision* all state vectors, i.e. particle coordinates, velocities and
forces, are stored in single precision, critical variables are double
precision. A typical example of the latter is the virial, which is a sum
over all forces in the system, which have varying signs. In addition, in
many parts of the code we managed to avoid double precision for
arithmetic, by paying attention to summation order or reorganization of
mathematical expressions. The default configuration uses mixed
precision, but it is easy to turn on double precision by adding the
option ``-DGMX\_DOUBLE=on`` to ``cmake``. Double
precision will be 20 to 100% slower than mixed precision depending on
the architecture you are running on. Double precision will use somewhat
more memory and run input, energy and full-precision trajectory files
will be almost twice as large.

The energies in mixed precision are accurate up to the last decimal, the
last one or two decimals of the forces are non-significant. The virial
is less accurate than the forces, since the virial is only one order of
magnitude larger than the size of each element in the sum over all atoms
(sec.
). In most cases this is not really a
problem, since the fluctuations in the virial can be two orders of
magnitude larger than the average. Using cut-offs for the Coulomb
interactions cause large errors in the energies, forces, and virial.
Even when using a reaction-field or lattice sum method, the errors are
larger than, or comparable to, the errors due to the partial use of
single precision. Since MD is chaotic, trajectories with very similar
starting conditions will diverge rapidly, the divergence is faster in
mixed precision than in double precision.

For most simulations, mixed precision is accurate enough. In some cases
double precision is required to get reasonable results:

-  normal mode analysis, for the conjugate gradient or l-bfgs
   minimization and the calculation and diagonalization of the Hessian

-  long-term energy conservation, especially for large systems

