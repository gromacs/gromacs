Molecule definition
-------------------

Moleculetype entries
~~~~~~~~~~~~~~~~~~~~

An organizational structure that usually corresponds to molecules is the
``[ moleculetype ]`` entry. This entry serves two main
purposes. One is to give structure to the topology file(s), usually
corresponding to real molecules. This makes the topology easier to read
and writing it less labor intensive. A second purpose is computational
efficiency. The system definition that is kept in memory is proportional
in size of the ``moleculetype`` definitions. If a molecule
is present in 100000 copies, this saves a factor of 100000 in memory,
which means the system usually fits in cache, which can improve
performance tremendously. Interactions that correspond to chemical
bonds, that generate exclusions, can only be defined between atoms
within a ``moleculetype``. It is allowed to have multiple
molecules which are not covalently bonded in one
``moleculetype`` definition. Molecules can be made
infinitely long by connecting to themselves over periodic boundaries.
When such periodic molecules are present, an option in the
:ref:`mdp` file needs to be set to tell |Gromacs| not to attempt
to make molecules that are broken over periodic boundaries whole again.

Intermolecular interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, one would like atoms in different molecules to also
interact with other interactions than the usual non-bonded interactions.
This is often the case in binding studies. When the molecules are
covalently bound, e.g. a ligand binding covalently to a protein, they
are effectively one molecule and they should be defined in one
``[ moleculetype ]`` entry. Note that
:ref:`pdb2gmx <gmx pdb2gmx>` has an option to put two or more molecules in
one ``[ moleculetype ]`` entry. When molecules are not
covalently bound, it is much more convenient to use separate
``moleculetype`` definitions and specify the intermolecular
interactions in the ``[ intermolecular_interactions]``
section. In this section, which is placed at the end of the topology
(see :numref:`Table %s <tab-topfile1>`), normal bonded interactions
can be specified using global atom indices. The only restrictions are
that no interactions can be used that generates exclusions and no
constraints can be used.

.. _pairinteractions:

Intramolecular pair interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extra Lennard-Jones and electrostatic interactions between pairs of
atoms in a molecule can be added in the ``[ pairs ]`` section of a molecule
definition. The parameters for these interactions can be set
independently from the non-bonded interaction parameters. In the GROMOS
force fields, pairs are only used to modify the 1-4 interactions
(interactions of atoms separated by three bonds). In these force fields
the 1-4 interactions are excluded from the non-bonded interactions (see
sec. :ref:`excl`).

::


    [ pairtypes ]
      ; i    j func         cs6          cs12 ; THESE ARE 1-4 INTERACTIONS
        O    O    1 0.22617E-02   0.74158E-06
        O   OM    1 0.22617E-02   0.74158E-06
        .....

The pair interaction parameters for the atom types in ``ffnonbonded.itp``
are listed in the ``[ pairtypes ]`` section. The GROMOS force fields list all these
interaction parameters explicitly, but this section might be empty for
force fields like OPLS that calculate the 1-4 interactions by uniformly
scaling the parameters. Pair parameters that are not present in the ``[ pairtypes ]``
section are only generated when ``gen-pairs`` is set to ``yes`` in the
``[ defaults ]`` directive of ``forcefield.itp`` (see :ref:`topfile`). When ``gen-pairs`` is
set to ``no``, :ref:`grompp <gmx grompp>` will give a warning for each pair type for which no
parameters are given.

The normal pair interactions, intended for 1-4 interactions, have
function type 1. Function type 2 and the ``[ pairs_nb ]`` are intended for free-energy
simulations. When determining hydration free energies, the solute needs
to be decoupled from the solvent. This can be done by adding a B-state
topology (see sec. :ref:`fecalc`) that uses zero for all solute
non-bonded parameters, *i.e.* charges and LJ parameters. However, the
free energy difference between the A and B states is not the total
hydration free energy. One has to add the free energy for reintroducing
the internal Coulomb and LJ interactions in the solute when in vacuum.
This second step can be combined with the first step when the Coulomb
and LJ interactions within the solute are not modified. For this
purpose, there is a pairs function type 2, which is identical to
function type 1, except that the B-state parameters are always identical
to the A-state parameters. For searching the parameters in the ``[ pairtypes ]`` section,
no distinction is made between function type 1 and 2. The pairs section
``[ pairs_nb ]`` is intended to replace the non-bonded interaction. It uses the unscaled
charges and the non-bonded LJ parameters; it also only uses the A-state
parameters. **Note** that one should add exclusions for all atom pairs
listed in ``[ pairs_nb ]``, otherwise such pairs will also end up in the normal neighbor
lists.

Alternatively, this same behavior can be achieved without ever touching
the topology, by using the ``couple-moltype``, ``couple-lambda0``,
``couple-lambda1``, and ``couple-intramol`` keywords. See sections
sec. :ref:`fecalc` and sec. :ref:`dgimplement` for more information.

All three pair types always use plain Coulomb interactions, even when
Reaction-field, PME, Ewald or shifted Coulomb interactions are selected
for the non-bonded interactions. Energies for types 1 and 2 are written
to the energy and log file in separate “LJ-14” and “Coulomb-14” entries
per energy group pair. Energies for ``[ pairs_nb ]`` are added to the “LJ-(SR)” and
“Coulomb-(SR)” terms.

.. _excl:

Exclusions
~~~~~~~~~~

The exclusions for non-bonded interactions are generated by :ref:`grompp <gmx grompp>` for
neighboring atoms up to a certain number of bonds away, as defined in
the ``[ moleculetype ]`` section in the topology file (see :ref:`topfile`). Particles are
considered bonded when they are connected by “chemical” bonds (``[ bonds ]`` types 1
to 5, 7 or 8) or constraints (``[ constraints ]`` type 1). Type 5 ``[ bonds ]`` can be used to create a
connection between two atoms without creating an interaction. There is a
harmonic interaction (``[ bonds ]`` type 6) that does not connect the atoms by a
chemical bond. There is also a second constraint type (``[ constraints ]`` type 2) that
fixes the distance, but does not connect the atoms by a chemical bond.
For a complete list of all these interactions, see :numref:`Table %s <tab-topfile2>`.

Extra exclusions within a molecule can be added manually in a
``[ exclusions ]`` section. Each line should start with one
atom index, followed by one or more atom indices. All non-bonded
interactions between the first atom and the other atoms will be
excluded.

When all non-bonded interactions within or between groups of atoms need
to be excluded, is it more convenient and much more efficient to use
energy monitor group exclusions (see sec. :ref:`groupconcept`).
