Shell molecular dynamics
------------------------

|Gromacs| can simulate polarizability using the shell model of Dick and
Overhauser \ :ref:`43 <refDick58>`. In such models a shell particle
representing the electronic degrees of freedom is attached to a nucleus
by a spring. The potential energy is minimized with respect to the shell
position at every step of the simulation (see below). Successful
applications of shell models in |Gromacs| have been published for
:math:`N_2` :ref:`44 <refJordan95>` and water\ :ref:`45 <refMaaren2001a>`.

Optimization of the shell positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The force :math:`\mathbf{F}`\ :math:`_S` on a shell
particle :math:`S` can be decomposed into two components

.. math:: \mathbf{F}_S ~=~ \mathbf{F}_{bond} + \mathbf{F}_{nb}
          :label: eqnshellforcedecomp

where :math:`\mathbf{F}_{bond}` denotes the
component representing the polarization energy, usually represented by a
harmonic potential and :math:`\mathbf{F}_{nb}` is the sum of Coulomb
and van der Waals interactions. If we assume that
:math:`\mathbf{F}_{nb}` is almost constant we
can analytically derive the optimal position of the shell, i.e. where
:math:`\mathbf{F}_S` = 0. If we have the shell S connected to atom A we have

.. math:: \mathbf{F}_{bond} ~=~ k_b \left( \mathbf{x}_S - \mathbf{x}_A\right).
          :label: eqnshell

In an iterative solver, we have positions :math:`\mathbf{x}_S(n)` where :math:`n` is
the iteration count. We now have at iteration :math:`n`

.. math:: \mathbf{F}_{nb} ~=~ \mathbf{F}_S - k_b \left( \mathbf{x}_S(n) - \mathbf{x}_A\right)
          :label: eqnshellsolv

and the optimal position for the shells :math:`x_S(n+1)` thus follows from

.. math:: \mathbf{F}_S - k_b \left( \mathbf{x}_S(n) - \mathbf{x}_A\right) + k_b \left( \mathbf{x}_S(n+1) - \mathbf{x}_A\right) = 0
          :label: eqnshelloptpos

if we write

.. math:: \Delta \mathbf{x}_S = \mathbf{x}_S(n+1) - \mathbf{x}_S(n)
          :label: eqnshelloptpos2

we finally obtain

.. math:: \Delta \mathbf{x}_S = \mathbf{F}_S/k_b
          :label: eqnshelloptpos3

which then yields the algorithm to compute the next trial in the
optimization of shell positions

.. math:: \mathbf{x}_S(n+1) ~=~ \mathbf{x}_S(n) + \mathbf{F}_S/k_b.
          :label: eqnshelloptpos4

