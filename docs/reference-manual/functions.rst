.. _ff:

Interaction function and force fields
=====================================

To accommodate the potential functions used in some popular force fields
(see :ref:`ff`), |Gromacs| offers a choice of functions, both for
non-bonded interaction and for dihedral interactions. They are described
in the appropriate subsections.

The potential functions can be subdivided into three parts

#. *Non-bonded*: Lennard-Jones or Buckingham, and Coulomb or modified
   Coulomb. The non-bonded interactions are computed on the basis of a
   neighbor list (a list of non-bonded atoms within a certain radius),
   in which exclusions are already removed.

#. *Bonded*: covalent bond-stretching, angle-bending, improper
   dihedrals, and proper dihedrals. These are computed on the basis of
   fixed lists.

#. *Restraints*: position restraints, angle restraints, distance
   restraints, orientation restraints and dihedral restraints, all based
   on fixed lists.

#. *Applied Forces*: externally applied forces, see
   chapter :ref:`special`.


.. toctree::
   :maxdepth: 2

   functions/nonbonded-interactions
   functions/bonded-interactions
   functions/restraints
   functions/polarization
   functions/free-energy-interactions
   functions/interaction-methods
   functions/long-range-electrostatics
   functions/long-range-vdw
   functions/force-field


