.. _defunits:

Definitions and Units
=====================

Notation
--------

The following conventions for mathematical typesetting are used
throughout this document:

.. |vecex| replace:: :math:`{\mathbf{r}_i}`
.. |lenex| replace:: :math:`r_i`

.. table:: 
    :align: center
    :widths: auto

    +---------------+-------------+---------+
    | Item          | Notation    | Example |
    +===============+=============+=========+
    | Vector        | Bold italic | |vecex| |
    +---------------+-------------+---------+
    | Vector Length | Italic      | |lenex| |
    +---------------+-------------+---------+

We define the *lowercase* subscripts :math:`i`, :math:`j`, :math:`k` and
:math:`l` to denote particles: :math:`\mathbf{r}_i` is the
*position vector* of particle :math:`i`, and using this notation:

.. math:: \begin{aligned}
          \mathbf{r}_{ij}	=	\mathbf{r}_j-\mathbf{r}_i\\
          r_{ij}=	| \mathbf{r}_{ij} | \end{aligned}
          :label: eqnnotation

The force on particle :math:`i` is denoted by
:math:`\mathbf{F}_i` and

.. math:: \mathbf{F}_{ij} = \mbox{force on $i$ exerted by $j$}
          :label: eqbforcenotation

MD units
--------

|Gromacs| uses a consistent set of units that produce values in the
vicinity of unity for most relevant molecular quantities. Let us call
them *MD units*. The basic units in this system are nm, ps, K, electron
charge (e) and atomic mass unit (u), see :numref:`Table %s <table-basicunits>`
The values used in |Gromacs| are
taken from the CODATA Internationally recommended 2010 values of
fundamental physical constants (see `NIST homepage <http://nist.gov>`__). 

.. |tnm| replace:: :math:`\mathrm{nm = }10^{-9}\ m`
.. |tu1| replace:: u (unified atomic mass unit) =
.. |tu2| replace:: :math:`1.660\,538\,921 \times 10^{-27}\ kg`
.. |tti| replace:: :math:`\mathrm{ps = }10^{-12}\ s`
.. |tc1| replace:: *e* = elementary charge =
.. |tc2| replace:: :math:`1.602\,176\,565 \times 10^{-19}\ C`
.. |tte| replace:: K 

.. _table-basicunits:

.. table:: Basic units used in |Gromacs|
    :align: center
    :widths: auto

    +--------------+--------+-------+
    | Quantity     | Symbol | Unit  |
    +==============+========+=======+
    | length       |     r  | |tnm| |
    +--------------+--------+-------+
    | mass         |     m  | |tu1| |
    |              |        | |tu2| |
    +--------------+--------+-------+
    | time         |     t  | |tti| |
    +--------------+--------+-------+
    | charge       |     q  | |tc1| |
    |              |        | |tc2| |
    +--------------+--------+-------+
    | temperature  |     T  | |tte| |
    +--------------+--------+-------+




Consistent
with these units are a set of derived units, given in
:numref:`Table %s <table-derivedunits>`

.. |tse|  replace:: :math:`E,V`
.. |tsf|  replace:: :math:`\mathbf{F}`
.. |tsp|  replace:: :math:`p`
.. |tsv|  replace:: :math:`v`
.. |tsd|  replace:: :math:`\mu`
.. |tsep| replace:: :math:`\Phi`
.. |tsef| replace:: :math:`E`
.. |tdue|   replace:: :math:`\mathrm{kJ~mol}^{-1}`
.. |tduf|   replace:: :math:`\mathrm{kJ~mol}^{-1}~\mathrm{nm}^{-1}`
.. |tdup|   replace:: bar
.. |tduv|   replace:: :math:`\mathrm{nm~ps}^{-1} = 1000\mathrm{~m~s}^{-1}`
.. |tdud|   replace:: :math:`\mathrm{e\ nm}`
.. |tduep1| replace:: :math:`\mathrm{kJ~mol}^{-1}\mathrm{~e}^{-1} =`
.. |tduep2| replace:: :math:`0.010\,364\,269\,19` Volt
.. |tduef1| replace:: :math:`\mathrm{kJ~mol}^{-1}\mathrm{~nm}^{-1}\ \mathrm{e}^{-1} =`
.. |tduef2| replace:: :math:`1.036\,426\,919 \times 10^7\mathrm{~V m}^{-1}`

.. _table-derivedunits:

.. table::
    Derived units. Note that an additional conversion factor of 10\ :math:`^{28}` a.m.u (\ :math:`\approx` 16.6)
    is applied to get bar instead of internal MD units in the energy and
    log files
    :align: center
    :widths: auto

    +--------------------+--------+----------+
    | Quantity           | Symbol | Unit     |
    +====================+========+==========+
    | energy             | |tse|  | |tdue|   |
    +--------------------+--------+----------+
    | Force              | |tsf|  | |tduf|   |
    +--------------------+--------+----------+
    | pressure           | |tsp|  | |tdup|   |
    +--------------------+--------+----------+
    | velocity           | |tsv|  | |tduv|   |
    +--------------------+--------+----------+
    | dipole moment      | |tsd|  | |tdud|   |
    +--------------------+--------+----------+
    | electric potential | |tsep| | |tduep1| |
    |                    |        | |tduep2| |
    +--------------------+--------+----------+
    | electric field     | |tsef| | |tduef1| |
    |                    |        | |tduef2| |
    +--------------------+--------+----------+


The **electric conversion factor**
:math:`f=\frac{1}{4 \pi \varepsilon_o}={138.935\,458}`
:math:`\mathrm{kJ}~\mathrm{mol}^{-1}\mathrm{nm}~\mathrm{ e}^{-2}`.
It relates the mechanical quantities to the electrical quantities as in

.. math:: V = f \frac{q^2}{r} \mbox{\ \ or\ \ } F = f \frac{q^2}{r^2}
          :label: eqnelecconv

Electric potentials :math:`\Phi` and electric fields
:math:`\mathbf{E}` are intermediate quantities in the
calculation of energies and forces. They do not occur inside |Gromacs|. If
they are used in evaluations, there is a choice of equations and related
units. We strongly recommend following the usual practice of including
the factor :math:`f` in expressions that evaluate :math:`\Phi` and
:math:`\mathbf{E}`:

.. math:: \begin{aligned}
          \Phi(\mathbf{r}) = f \sum_j \frac{q_j}{| \mathbf{r}-\mathbf{r}_j | } 	\\
          \mathbf{E}(\mathbf{r}) = f \sum_j q_j \frac{(\mathbf{r}-\mathbf{r}_j)}{| \mathbf{r}-\mathbf{r}_j| ^3}\end{aligned}
          :label: eqnelecfacinclude

With these definitions, :math:`q\Phi` is an energy and
:math:`q\mathbf{E}` is a force. The units are those given
in :numref:`Table %s <table-derivedunits>`
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
values (see :numref:`Table %s <table-consts>`). All quantities are per
mol rather than per molecule. There is no distinction between
Boltzmann’s constant :math:`k` and the gas constant :math:`R`: their
value is :math:`0.008\,314\,462\,1\mathrm{kJ~mol}^{-1} \mathrm{K}^{-1}`.

.. _table-consts:

.. table:: 
    Some Physical Constants
    :align: center
    :widths: auto

    +----------------+----------------------+--------------------------------------------------------------------------+
    | Symbol         | Name                 | Value                                                                    |
    +================+======================+==========================================================================+
    | :math:`N_{AV}` | Avogadro's number    | :math:`6.022\,141\,29\times 10^{23}~\mathrm{mol}^{-1}`                   |
    +----------------+----------------------+--------------------------------------------------------------------------+
    | :math:`R`      | gas constant         | :math:`8.314\,462\,1\times 10^{-3}~\mathrm{kJ~mol}^{-1}~\mathrm{K}^{-1}` |
    +----------------+----------------------+--------------------------------------------------------------------------+
    | :math:`k_B`    | Boltzmann's constant | *idem*                                                                   |
    +----------------+----------------------+--------------------------------------------------------------------------+
    | :math:`h`      | Planck's constant    | :math:`0.399\,031\,271~\mathrm{kJ~mol}^{-1}~\mathrm{ps}`                 |
    +----------------+----------------------+--------------------------------------------------------------------------+
    | :math:`\hbar`  | Dirac's constant     | :math:`0.063\,507\,799\,3~\mathrm{kJ~mol}^{-1}~\mathrm{ps}`              |
    +----------------+----------------------+--------------------------------------------------------------------------+
    | :math:`c`      | velocity of light    | :math:`299\,792.458~\mathrm{nm~ps}^{-1}`                                 |
    +----------------+----------------------+--------------------------------------------------------------------------+



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
reduced temperature. A |Gromacs| temperature :math:`T=1` means a reduced
temperature of :math:`0.008\ldots` units; if a reduced temperature of 1
is required, the |Gromacs| temperature should be :math:`120.272\,36`.

In :numref:`Table %s <table-reduced>` quantities are given for LJ
potentials:

.. math:: V_{LJ} = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]
          :label: eqnbaseljpotentials

.. _table-reduced:

.. table:: 
    Reduced Lennard-Jones quantities
    :align: center
    :widths: auto

    +-------------+----------------+------------------------------------------+
    | Quantity    | Symbol         | Relation to SI                           |
    +=============+================+==========================================+
    | Length      | r\ :math:`^*`  | r\ :math:`\sigma^{-1}`                   |
    +-------------+----------------+------------------------------------------+
    | Mass        | m\ :math:`^*`  | m M\ :math:`^{-1}`                       |
    +-------------+----------------+------------------------------------------+
    | Time        | t\ :math:`^*`  | t\ :math:`\sigma^{-1}~\sqrt{\epsilon/M}` |
    +-------------+----------------+------------------------------------------+
    | Temperature | T\ :math:`^*`  | k\ :math:`_B\mathrm{T}~\epsilon^{-1}`    |
    +-------------+----------------+------------------------------------------+
    | Energy      | E\ :math:`^*`  | E\ :math:`\epsilon^{-1}`                 |
    +-------------+----------------+------------------------------------------+
    | Force       | F\ :math:`^*`  | F\ :math:`\sigma~\epsilon^{-1}`          |
    +-------------+----------------+------------------------------------------+
    | Pressure    | P\ :math:`^*`  | P\ :math:`\sigma ^3 \epsilon^{-1}`       |
    +-------------+----------------+------------------------------------------+
    | Velocity    | v\ :math:`^*`  | v\ :math:`\sqrt{M/\epsilon}`             |
    +-------------+----------------+------------------------------------------+
    | Density     | :math:`\rho^*` | N\ :math:`\sigma ^3~V^{-1}`              |
    +-------------+----------------+------------------------------------------+




Mixed or Double precision
-------------------------

|Gromacs| can be compiled in either mixed or double precision.
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
option ``-DGMX\_DOUBLE=on`` to ``cmake``. Double precision will be 20 to 100%
slower than mixed precision depending on the architecture you are
running on. Double precision will use somewhat more memory and run
input, energy and full-precision trajectory files will be almost twice
as large.

The energies in mixed precision are accurate up to the last decimal, the
last one or two decimals of the forces are non-significant. The virial
is less accurate than the forces, since the virial is only one order of
magnitude larger than the size of each element in the sum over all atoms
(sec. :ref:`virial`). In most cases this is not really a problem, since
the fluctuations in the virial can be two orders of magnitude larger
than the average. Using cut-offs for the Coulomb interactions cause
large errors in the energies, forces, and virial. Even when using a
reaction-field or lattice sum method, the errors are larger than, or
comparable to, the errors due to the partial use of single precision.
Since MD is chaotic, trajectories with very similar starting conditions
will diverge rapidly, the divergence is faster in mixed precision than
in double precision.

For most simulations, mixed precision is accurate enough. In some cases
double precision is required to get reasonable results:

-  normal mode analysis, for the conjugate gradient or l-bfgs
   minimization and the calculation and diagonalization of the Hessian

-  long-term energy conservation, especially for large systems

.. raw:: latex

    \clearpage


