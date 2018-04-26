.. _fepmf:

Calculating a PMF using the free-energy code
--------------------------------------------

The free-energy coupling-parameter approach (see sec. :ref:`fecalc`)
provides several ways to calculate potentials of mean force. A potential
of mean force between two atoms can be calculated by connecting them
with a harmonic potential or a constraint. For this purpose there are
special potentials that avoid the generation of extra exclusions,
see sec. :ref:`excl`. When the position of the minimum or the constraint
length is 1 nm more in state B than in state A, the restraint or
constraint force is given by :math:`\partial H/\partial \lambda`. The
distance between the atoms can be changed as a function of
:math:`\lambda` and time by setting delta-lambda in the :ref:`mdp` file. The
results should be identical (although not numerically due to the
different implementations) to the results of the pull code with umbrella
sampling and constraint pulling. Unlike the pull code, the free energy
code can also handle atoms that are connected by constraints.

Potentials of mean force can also be calculated using position
restraints. With position restraints, atoms can be linked to a position
in space with a harmonic potential (see :ref:`positionrestraint`).
These positions can be made a function of the coupling parameter
:math:`\lambda`. The positions for the A and the B states are supplied
to :ref:`grompp <gmx grompp>` with the ``-r`` and ``-rb`` options, respectively. One could use this
approach to do targeted MD; note that we do not encourage the use of
targeted MD for proteins. A protein can be forced from one conformation
to another by using these conformations as position restraint
coordinates for state A and B. One can then slowly change
:math:`\lambda` from 0 to 1. The main drawback of this approach is that
the conformational freedom of the protein is severely limited by the
position restraints, independent of the change from state A to B. Also,
the protein is forced from state A to B in an almost straight line,
whereas the real pathway might be very different. An example of a more
fruitful application is a solid system or a liquid confined between
walls where one wants to measure the force required to change the
separation between the boundaries or walls. Because the boundaries (or
walls) already need to be fixed, the position restraints do not limit
the system in its sampling.
