.. _groupconcept:

The group concept
-----------------

The |Gromacs| MD and analysis programs use user-defined *groups* of atoms
to perform certain actions on. The maximum number of groups is 256, but
each atom can only belong to six different groups, one each of the
following:

temperature-coupling group
    The temperature coupling parameters (reference temperature, time
    constant, number of degrees of freedom, see :ref:`update`) can be
    defined for each T-coupling group separately. For example, in a
    solvated macromolecule the solvent (that tends to generate more
    heating by force and integration errors) can be coupled with a
    shorter time constant to a bath than is a macromolecule, or a
    surface can be kept cooler than an adsorbing molecule. Many
    different T-coupling groups may be defined. See also center of mass
    groups below.

freeze group

    Atoms that belong to a freeze group are kept stationary in the
    dynamics. This is useful during equilibration, *e.g.* to avoid badly
    placed solvent molecules giving unreasonable kicks to protein atoms,
    although the same effect can also be obtained by putting a
    restraining potential on the atoms that must be protected. The
    freeze option can be used, if desired, on just one or two
    coordinates of an atom, thereby freezing the atoms in a plane or on
    a line. When an atom is partially frozen, constraints will still be
    able to move it, even in a frozen direction. A fully frozen atom can
    not be moved by constraints. Many freeze groups can be defined.
    Frozen coordinates are unaffected by pressure scaling; in some cases
    this can produce unwanted results, particularly when constraints are
    also used (in this case you will get very large pressures).
    Accordingly, it is recommended to avoid combining freeze groups with
    constraints and pressure coupling. For the sake of equilibration it
    could suffice to start with freezing in a constant volume
    simulation, and afterward use position restraints in conjunction
    with constant pressure.

accelerate group

    On each atom in an “accelerate group” an acceleration
    :math:`\mathbf{a}^g` is imposed. This is equivalent to
    an external force. This feature makes it possible to drive the
    system into a non-equilibrium state and enables the performance of
    non-equilibrium MD and hence to obtain transport properties.

energy-monitor group

    Mutual interactions between all energy-monitor groups are compiled
    during the simulation. This is done separately for Lennard-Jones and
    Coulomb terms. In principle up to 256 groups could be defined, but
    that would lead to 256\ :math:`\times`\ 256 items! Better use this
    concept sparingly.

    All non-bonded interactions between pairs of energy-monitor groups
    can be excluded (see details in the User Guide). Pairs of particles
    from excluded pairs of energy-monitor groups are not put into the
    pair list. This can result in a significant speedup for simulations
    where interactions within or between parts of the system are not
    required.

center of mass group

    In |Gromacs|, the center of mass (COM) motion can be removed, for
    either the complete system or for groups of atoms. The latter is
    useful, *e.g.* for systems where there is limited friction (*e.g.*
    gas systems) to prevent center of mass motion to occur. It makes
    sense to use the same groups for temperature coupling and center of
    mass motion removal.

Compressed position output group

    In order to further reduce the size of the compressed trajectory
    file (:ref:`xtc` or :ref:`tng`), it is possible to
    store only a subset of all particles. All x-compression groups that
    are specified are saved, the rest are not. If no such groups are
    specified, than all atoms are saved to the compressed trajectory
    file.

The use of groups in |Gromacs| tools is described in
sec. :ref:`usinggroups`.
