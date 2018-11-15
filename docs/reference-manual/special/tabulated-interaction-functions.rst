Tabulated interaction functions
-------------------------------

.. _cubicspline:

Cubic splines for potentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some of the inner loops of |Gromacs|, look-up tables are used for
computation of potential and forces. The tables are interpolated using a
cubic spline algorithm. There are separate tables for electrostatic,
dispersion, and repulsion interactions, but for the sake of caching
performance these have been combined into a single array. The cubic
spline interpolation for :math:`x_i \leq x < x_{i+1}` looks like this:

.. math::  V_s(x) = A_0 + A_1 \,\epsilon + A_2 \,\epsilon^2 + A_3 \,\epsilon^3
           :label: eqnspline

where the table spacing :math:`h` and fraction :math:`\epsilon` are
given by:

.. math::  \begin{aligned}
           h	&=&	x_{i+1} - x_i	\\
           \epsilon&=&	(x - x_i)/h\end{aligned}
           :label: eqntablespaceing

so that :math:`0 \le \epsilon < 1`. From this, we can calculate the
derivative in order to determine the forces:

.. math::  -V_s'(x) ~=~ 
           -\frac{{\rm d}V_s(x)}{{\rm d}\epsilon}\frac{{\rm d}\epsilon}{{\rm d}x} ~=~
           -(A_1 + 2 A_2 \,\epsilon + 3 A_3 \,\epsilon^2)/h
           :label: eqntablederivative

The four coefficients are determined from the four conditions that
:math:`V_s` and :math:`-V_s'` at both ends of each interval should match
the exact potential :math:`V` and force :math:`-V'`. This results in the
following errors for each interval:

.. math:: \begin{aligned}
          | V_s  - V  | _{max} &=& V'''' \frac{h^4}{384} + O(h^5) \\
          | V_s' - V' | _{max} &=& V'''' \frac{h^3}{72\sqrt{3}} + O(h^4) \\
          | V_s''- V''| _{max} &=& V'''' \frac{h^2}{12}  + O(h^3)\end{aligned}
          :label: eqntableerrors

V and V’ are continuous, while V” is the first discontinuous
derivative. The number of points per nanometer is 500 and 2000 for
mixed- and double-precision versions of |Gromacs|, respectively. This
means that the errors in the potential and force will usually be smaller
than the mixed precision accuracy.

|Gromacs| stores :math:`A_0`, :math:`A_1`, :math:`A_2` and :math:`A_3`.
The force routines get a table with these four parameters and a scaling
factor :math:`s` that is equal to the number of points per nm. (**Note**
that :math:`h` is :math:`s^{-1}`). The algorithm goes a little something
like this:

#. Calculate distance vector
   (:math:`\mathbf{r}_{ij}`) and distance
   r\ :math:`_{ij}`

#. Multiply r\ :math:`_{ij}` by :math:`s` and truncate to an integer
   value :math:`n_0` to get a table index

#. Calculate fractional component (:math:`\epsilon` =
   :math:`s`\ r\ :math:`_{ij} - n_0`) and :math:`\epsilon^2`

#. Do the interpolation to calculate the potential :math:`V` and the
   scalar force :math:`f`

#. Calculate the vector force :math:`\mathbf{F}` by
   multiplying :math:`f` with
   :math:`\mathbf{r}_{ij}`

**Note** that table look-up is significantly *slower* than computation
of the most simple Lennard-Jones and Coulomb interaction. However, it is
much faster than the shifted Coulomb function used in conjunction with
the PPPM method. Finally, it is much easier to modify a table for the
potential (and get a graphical representation of it) than to modify the
inner loops of the MD program.

User-specified potential functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use your own potential functions without editing the
|Gromacs| code. The potential function should be according to the
following equation

.. math:: V(r_{ij}) ~=~ \frac{q_i q_j}{4 \pi\epsilon_0} f(r_{ij}) + C_6 \,g(r_{ij}) + C_{12} \,h(r_{ij})
          :label: eqnuserpotfunction

where :math:`f`, :math:`g`, and :math:`h` are user defined functions.
**Note** that if :math:`g(r)` represents a normal dispersion
interaction, :math:`g(r)` should be :math:`<` 0. C\ :math:`_6`,
C\ :math:`_{12}` and the charges are read from the topology. Also note
that combination rules are only supported for Lennard-Jones and
Buckingham, and that your tables should match the parameters in the
binary topology.

When you add the following lines in your :ref:`mdp` file:

::

    rlist           = 1.0
    coulombtype     = User
    rcoulomb        = 1.0
    vdwtype         = User
    rvdw            = 1.0

:ref:`mdrun <gmx mdrun>` will read a single non-bonded table file, or
multiple when ``energygrp-table`` is set (see below). The
name of the file(s) can be set with the :ref:`mdrun <gmx mdrun>` option
``-table``. The table file should contain seven columns of
table look-up data in the order: :math:`x`, :math:`f(x)`,
:math:`-f'(x)`, :math:`g(x)`, :math:`-g'(x)`, :math:`h(x)`,
:math:`-h'(x)`. The :math:`x` should run from 0 to :math:`r_c+1` (the
value of ``table_extension`` can be changed in the :ref:`mdp` file). You can
choose the spacing you like; for the standard tables |Gromacs| uses a
spacing of 0.002 and 0.0005 nm when you run in mixed and double
precision, respectively. In this context, :math:`r_c` denotes the
maximum of the two cut-offs ``rvdw`` and ``rcoulomb`` (see above). These
variables need not be the same (and need not be 1.0 either). Some
functions used for potentials contain a singularity at :math:`x = 0`,
but since atoms are normally not closer to each other than 0.1 nm, the
function value at :math:`x = 0` is not important. Finally, it is also
possible to combine a standard Coulomb with a modified LJ potential (or
vice versa). One then specifies *e.g.* ``coulombtype = Cut-off`` or
``coulombtype = PME``, combined with ``vdwtype = User``. The table file must
always contain the 7 columns however, and meaningful data (i.e. not
zeroes) must be entered in all columns. A number of pre-built table
files can be found in the ``GMXLIB`` directory for 6-8, 6-9, 6-10, 6-11, and
6-12 Lennard-Jones potentials combined with a normal Coulomb.

If you want to have different functional forms between different groups
of atoms, this can be set through energy groups. Different tables can be
used for non-bonded interactions between different energy groups pairs
through the :ref:`mdp` option ``energygrp-table`` (see details in the User Guide).
Atoms that should interact with a different potential should be put into
different energy groups. Between group pairs which are not listed in
``energygrp-table``, the normal user tables will be used. This makes it easy
to use a different functional form between a few types of atoms.
