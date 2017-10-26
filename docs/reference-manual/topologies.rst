Topologies
==========

Introduction
------------

|Gromacs| must know on which atoms and combinations of atoms the various
contributions to the potential functions (see chapter :ref:`ff`) must act.
It must also know what parameters must be applied to the various
functions. All this is described in the *topology* file :ref:`top`, which
lists the *constant attributes* of each atom. There are many more atom
types than elements, but only atom types present in biological systems
are parameterized in the force field, plus some metals, ions and
silicon. The bonded and special interactions are determined by fixed
lists that are included in the topology file. Certain non-bonded
interactions must be excluded (first and second neighbors), as these are
already treated in bonded interactions. In addition, there are *dynamic
attributes* of atoms - their positions, velocities and forces. These do
not strictly belong to the molecular topology, and are stored in the
coordinate file :ref:`gro` (positions and velocities), or
trajectory file :ref:`trr` (positions, velocities, forces).

This chapter describes the setup of the topology file, the :ref:`top` file and
the database files: what the parameters stand for and how/where to
change them if needed. First, all file formats are explained. Section
:ref:`fffiles` describes the organization of the files in each force
field.

**Note:** if you construct your own topologies, we encourage you to
upload them to our topology archive at our `webpage`_! Just imagine how thankful
you’d have been if your topology had been available there before you
started. The same goes for new force fields or modified versions of the
standard force fields - contribute them to the force field archive!

.. _homepage: `webpage`_

Particle type
-------------

In |Gromacs|, there are three types of
particles
, see :numref:`Table %s <tab-ptype>`. Only regular atoms and virtual
interaction sites are used in |Gromacs|; shells are necessary for
polarizable models like the Shell-Water models \ :ref:`45 <refMaaren2001a>`.

.. _tab-ptype:

.. table:: Particle types in |Gromacs|

           +--------------+----------+
           | Particle     | Symbol   |
           +==============+==========+
           | atom         | A        |
           +--------------+----------+
           | shell        | S        |
           +--------------+----------+
           | virtual side | V (or D) |
           +--------------+----------+


.. _atomtype:

Atom types
~~~~~~~~~~

Each force field defines a set of atom
types,
which have a characteristic name or number, and mass (in a.m.u.). These
listings are found in the ``atomtypes.atp`` file (:ref:`atp` =
**a**\ tom **t**\ ype **p**\ arameter file). Therefore, it is in this
file that you can begin to change and/or add an atom type. A sample from
the ``gromos43a1.ff`` force field is listed below.

::

     |  O  15.99940 ;     carbonyl oxygen (C=O)
     | OM  15.99940 ;     carboxyl oxygen (CO-)
     | OA  15.99940 ;     hydroxyl, sugar or ester oxygen
     | OW  15.99940 ;     water oxygen
     |  N  14.00670 ;     peptide nitrogen (N or NH)
     | NT  14.00670 ;     terminal nitrogen (NH2)
     | NL  14.00670 ;     terminal nitrogen (NH3)
     | NR  14.00670 ;     aromatic nitrogen
     | NZ  14.00670 ;     Arg NH (NH2)
     | NE  14.00670 ;     Arg NE (NH)
     |  C  12.01100 ;     bare carbon
     |CH1  13.01900 ;     aliphatic or sugar CH-group
     |CH2  14.02700 ;     aliphatic or sugar CH2-group
     |CH3  15.03500 ;     aliphatic CH3-group

**Note:** |Gromacs| makes use of the atom types as a name, *not* as a
number (as *e.g.* in GROMOS).

.. _vsitetop:

Virtual sites
~~~~~~~~~~~~~

Some force fields use virtual interaction sites (interaction sites that
are constructed from other particle positions) on which certain
interactions are located (*e.g.* on benzene rings, to reproduce the
correct quadrupole). This is described in sec. :ref:`virtualsites`.

To make virtual sites in your system, you should include a section
``[ virtual_sites? ]`` (for backward compatibility the old
name ``[ dummies? ]`` can also be used) in your topology
file, where the ``?`` stands for the number constructing
particles for the virtual site. This will be `:ref:`2`` for
type 2, `:ref:`3`` for types 3, 3fd, 3fad and 3out and
`:ref:`4`` for type 4fdn. The last of these replace an older
4fd type (with the ‘type’ value 1) that could occasionally be unstable;
while it is still supported internally in the code, the old 4fd type
should not be used in new input files. The different types are explained
in sec. :ref:`virtualsites`.

Parameters for type 2 should look like this:

::

    [ virtual_sites2 ]
    ; Site  from        funct  a
    5       1     2     1      0.7439756

for type 3 like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          b
    5       1     2     3      1       0.7439756  0.128012

for type 3fd like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          d
    5       1     2     3      2       0.5        -0.105

for type 3fad like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   theta      d
    5       1     2     3      3       120        0.5

for type 3out like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          b          c
    5       1     2     3      4       -0.4       -0.4       6.9281

for type 4fdn like this:

::

    [ virtual_sites4 ]
    ; Site  from                      funct   a          b          c
    5       1     2     3     4       2       1.0        0.9       0.105

This will result in the construction of a virtual site, number 5 (first
column ``Site``), based on the positions of the atoms
whose indices are 1 and 2 or 1, 2 and 3 or 1, 2, 3 and 4 (next two,
three or four columns ``from``) following the rules
determined by the function number (next column ``funct``)
with the parameters specified (last one, two or three columns
``a b . .``). Obviously, the atom numbers (including
virtual site number) depend on the molecule. It may be instructive to
study the topologies for TIP4P or TIP5P water models that are included
with the |Gromacs| distribution.

**Note** that if any constant bonded interactions are defined between
virtual sites and/or normal atoms, they will be removed by
:ref:`grompp <gmx grompp>` (unless the option ``-normvsbds`` is used). This
removal of bonded interactions is done after generating exclusions, as
the generation of exclusions is based on “chemically” bonded
interactions.

Virtual sites can be constructed in a more generic way using basic
geometric parameters. The directive that can be used is ``[ virtual_sitesn ]``. Required
parameters are listed in :numref:`Table %s <tab-topfile2>`. An example entry for
defining a virtual site at the center of geometry of a given set of
atoms might be:

::

    [ virtual_sitesn ]
    ; Site   funct    from
    5        1        1     2     3     4

Parameter files
---------------

Atoms
~~~~~

The *static* properties (see  :numref:`Table %s <tab-statprop>`) assigned to the atom
types are assigned based on data in several places. The mass is listed
in ``atomtypes.atp`` (see :ref:`atomtype`), whereas the charge is listed
in :ref:`rtp` (:ref:`rtp` = **r**\ esidue **t**\ opology **p**\ arameter file,
see :ref:`rtp`). This implies that the charges are only defined in the
building blocks of amino acids, nucleic acids or otherwise, as defined
by the user. When generating a :ref:`topology <top>` using the :ref:`pdb2gmx <gmx pdb2gmx>`
program, the information from these files is combined.

.. _tab-statprop:

.. table:: Static atom type properties in |Gromacs|

           +----------+------------------+----------+
           | Property | Symbol           | Unit     |
           +==========+==================+==========+
           | Type     | -                | -        |
           +----------+------------------+----------+
           | Mass     | m                | a.m.u.   |
           +----------+------------------+----------+
           | Charge   | q                | electron |
           +----------+------------------+----------+
           | epsilon  | :math:`\epsilon` | kJ/mol   |
           +----------+------------------+----------+
           | sigma    | :math:`\sigma`   | nm       |
           +----------+------------------+----------+


.. _nbpar:

Non-bonded parameters
~~~~~~~~~~~~~~~~~~~~~

The non-bonded parameters consist of the van der Waals parameters V (``c6``
or :math:`\sigma`, depending on the combination rule) and W (``c12`` or
:math:`\epsilon`), as listed in the file ``ffnonbonded.itp``, where ``ptype`` is
the particle type (see :numref:`Table %s <tab-ptype>`). As with the bonded
parameters, entries in ``[ *type ]`` directives are applied to their counterparts in
the topology file. Missing parameters generate warnings, except as noted
below in section :ref:`pairinteractions`.

::

    [ atomtypes ]
    ;name   at.num      mass      charge   ptype         V(c6)        W(c12)
        O        8  15.99940       0.000       A   0.22617E-02   0.74158E-06
       OM        8  15.99940       0.000       A   0.22617E-02   0.74158E-06
       .....

    [ nonbond_params ]
      ; i    j func       V(c6)        W(c12)
        O    O    1 0.22617E-02   0.74158E-06
        O   OA    1 0.22617E-02   0.13807E-05
        .....

**Note** that most of the included force fields also include the ``at.num.``
column, but this same information is implied in the OPLS-AA ``bond_type``
column. The interpretation of the parameters V and W depends on the
combination rule that was chosen in the ``[ defaults ]`` section of the topology file
(see :ref:`topfile`):

.. math::

   \begin{aligned}
   \mbox{for combination rule 1}: & &
   \begin{array}{llllll}
     \mbox{V}_{ii} & = & C^{(6)}_{i}  & = & 4\,\epsilon_i\sigma_i^{6} &
     \mbox{[ kJ mol$^{-1}$ nm$^{6}$ ]}\\
     \mbox{W}_{ii} & = & C^{(12)}_{i} & = & 4\,\epsilon_i\sigma_i^{12} &
     \mbox{[ kJ mol$^{-1}$ nm$^{12}$ ]}\\
   \end{array}
   \\
   \mbox{for combination rules 2 and 3}: & &
   \begin{array}{llll}
     \mbox{V}_{ii} & = & \sigma_i   & \mbox{[ nm ]} \\
     \mbox{W}_{ii} & = & \epsilon_i & \mbox{[ kJ mol$^{-1}$ ]}
   \end{array}\end{aligned}

Some or all combinations for different atom types can be given in the
``[ nonbond_params ]`` section, again with parameters V and
W as defined above. Any combination that is not given will be computed
from the parameters for the corresponding atom types, according to the
combination rule:

.. math::

   \begin{aligned}
   \mbox{for combination rules 1 and 3}: & &
   \begin{array}{lll}
     C^{(6)}_{ij}  & = & \left(C^{(6)}_i\,C^{(6)}_j\right)^{\frac{1}{2}} \\
     C^{(12)}_{ij} & = & \left(C^{(12)}_i\,C^{(12)}_j\right)^{\frac{1}{2}}
   \end{array}
   \\
   \mbox{for combination rule 2}: & &
   \begin{array}{lll}
     \sigma_{ij}   & = & \frac{1}{2}(\sigma_i+\sigma_j) \\
     \epsilon_{ij} & = & \sqrt{\epsilon_i\,\epsilon_j}
   \end{array}\end{aligned}

When :math:`\sigma` and :math:`\epsilon` need to be supplied (rules 2
and 3), it would seem it is impossible to have a non-zero :math:`C^{12}`
combined with a zero :math:`C^6` parameter. However, providing a
negative :math:`\sigma` will do exactly that, such that :math:`C^6` is
set to zero and :math:`C^{12}` is calculated normally. This situation
represents a special case in reading the value of :math:`\sigma`, and
nothing more.

There is only one set of combination rules for Buckingham potentials:

.. math::

   \begin{array}{rcl}
   A_{ij}   &=& \left(A_{ii} \, A_{jj}\right)^{1/2}    \\
   B_{ij}   &=& 2 / \left(\frac{1}{B_{ii}} + \frac{1}{B_{jj}}\right)        \\
   C_{ij}   &=& \left(C_{ii} \, C_{jj}\right)^{1/2}
   \end{array}

Bonded parameters
~~~~~~~~~~~~~~~~~

The bonded
parameters
(*i.e.* bonds, bond angles, improper and proper dihedrals) are listed in
``ffbonded.itp``.  The entries in this database describe,
respectively, the atom types in the interactions, the type of the
interaction, and the parameters associated with that interaction. These
parameters are then read by
:ref:`grompp <gmx grompp>` when processing a
topology and applied to the relevant bonded parameters, *i.e.*
``bondtypes`` are applied to entries in the
``[ bonds ]`` directive, etc. Any bonded parameter that is
missing from the relevant :``[ *type ]`` directive generates
a fatal error. The types of interactions are listed in
:numref:`Table %s <tab-topfile2>`. Example excerpts from such files
follow:

::

    [ bondtypes ]
      ; i    j func        b0          kb
        C    O    1   0.12300     502080.
        C   OM    1   0.12500     418400.
        ......

    [ angletypes ]
      ; i    j    k func       th0         cth
       HO   OA    C    1   109.500     397.480
       HO   OA  CH1    1   109.500     397.480
       ......

    [ dihedraltypes ]
      ; i    l func        q0          cq
     NR5*  NR5    2     0.000     167.360
     NR5* NR5*    2     0.000     167.360
     ......

    [ dihedraltypes ]
      ; j    k func      phi0          cp   mult
        C   OA    1   180.000      16.736      2
        C    N    1   180.000      33.472      2
        ......

    [ dihedraltypes ]
    ;
    ; Ryckaert-Bellemans Dihedrals
    ;
    ; aj    ak      funct
    CP2     CP2     3       9.2789  12.156  -13.120 -3.0597 26.240  -31.495

In the ``ffbonded.itp`` file, you can add bonded parameters.
If you want to include parameters for new atom types, make sure you
define them in ``atomtypes.atp`` as well.

For most interaction types, bonded parameters are searched and assigned
using an exact match for all type names and allowing only a single set
of parameters. The exception to this rule are
dihedral
parameters. For
``[ dihedraltypes ]`` wildcard atom type names can be
specified with the letter ``X`` in one or more of the four
positions. Thus one can for example assign proper dihedral parameters
based on the types of the middle two atoms. The parameters for the entry
with the most exact matches, i.e. the least wildcard matches, will be
used. Note that |Gromacs| versions older than 5.1.3 used the first match,
which means that a full match would be ignored if it is preceded by an
entry that matches on wildcards. Thus it is suggested to put wildcard
entries at the end, in case someone might use a forcefield with older
versions of |Gromacs|. In addition there is a dihedral type 9 which adds
the possibility of assigning multiple dihedral potentials, useful for
combining terms with different multiplicities. The different dihedral
potential parameter sets should be on directly adjacent lines in the
``[ dihedraltypes ]`` section.

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

.. _constraintalg:

Constraint algorithms
---------------------

Constraints are defined in the ``[ constraints ]`` section. The format is two atom numbers
followed by the function type, which can be 1 or 2, and the constraint
distance. The only difference between the two types is that type 1 is
used for generating exclusions and type 2 is not (see sec. :ref:`excl`).
The distances are constrained using the LINCS or the SHAKE algorithm,
which can be selected in the :ref:`mdp` file. Both types of constraints can be
perturbed in free-energy calculations by adding a second constraint
distance (see :ref:`constraintforce`). Several types of bonds and
angles (see :numref:`Table %s <tab-topfile2>`) can be converted automatically to
constraints by :ref:`grompp <gmx grompp>`. There are several options for this in the :ref:`mdp`
file.

We have also implemented the SETTLE
algorithm \ :ref:`47 <refMiyamoto92>`, which is an analytical solution of SHAKE, specifically for
water. SETTLE can be selected in the topology file. See, for instance,
the SPC molecule definition:

::

    [ moleculetype ]
    ; molname       nrexcl
    SOL             1

    [ atoms ]
    ; nr    at type res nr  ren nm  at nm   cg nr   charge
    1       OW      1       SOL     OW1     1       -0.82
    2       HW      1       SOL     HW2     1        0.41
    3       HW      1       SOL     HW3     1        0.41

    [ settles ]
    ; OW    funct   doh     dhh
    1       1       0.1     0.16333

    [ exclusions ]
    1       2       3
    2       1       3
    3       1       2

The ``[ settles ]`` directive defines the first atom of the
water molecule. The settle funct is always 1, and the distance between
O-H and H-H distances must be given. **Note** that the algorithm can
also be used for TIP3P and TIP4P \ :ref:`128 <refJorgensen83>`. TIP3P just has
another geometry. TIP4P has a virtual site, but since that is generated
it does not need to be shaken (nor stirred).

.. _pdb2gmxfiles:

:ref:`pdb2gmx <gmx pdb2gmx>` input files
----------------------------------------

The |Gromacs| program :ref:`pdb2gmx <gmx pdb2gmx>` generates a topology for the input
coordinate file. Several formats are supported for that coordinate file,
but :ref:`pdb` is the most commonly-used format (hence the name :ref:`pdb2gmx <gmx pdb2gmx>`).
:ref:`pdb2gmx <gmx pdb2gmx>` searches for force fields in sub-directories of the |Gromacs|
``share/top`` directory and your working directory. Force fields are
recognized from the file ``forcefield.itp`` in a directory with the
extension ``.ff``. The file ``forcefield.doc`` may be present, and if so, its
first line will be used by :ref:`pdb2gmx <gmx pdb2gmx>` to present a short description to the
user to help in choosing a force field. Otherwise, the user can choose a
force field with the ``-ff xxx`` command-line argument to :ref:`pdb2gmx <gmx pdb2gmx>`, which
indicates that a force field in a ``xxx.ff`` directory is desired. :ref:`pdb2gmx <gmx pdb2gmx>`
will search first in the working directory, then in the |Gromacs|
``share/top`` directory, and use the first matching ``xxx.ff`` directory found.

Two general files are read by :ref:`pdb2gmx <gmx pdb2gmx>`: an atom type file (extension
:ref:`atp`, see :ref:`atomtype`) from the force-field directory, and a file
called ``residuetypes.dat`` from either the working directory, or the
|Gromacs| ``share/top`` directory. ``residuetypes.dat`` determines which residue
names are considered protein, DNA, RNA, water, and ions.

:ref:`pdb2gmx <gmx pdb2gmx>` can read one or multiple databases with topological information
for different types of molecules. A set of files belonging to one
database should have the same basename, preferably telling something
about the type of molecules (*e.g.* aminoacids, rna, dna). The possible
files are:

-  ``<basename>.rtp``

-  ``<basename>.r2b (optional)``

-  ``<basename>.arn (optional)``

-  ``<basename>.hdb (optional)``

-  ``<basename>.n.tdb (optional)``

-  ``<basename>.c.tdb (optional)``

Only the :ref:`rtp` file, which contains the topologies of the building
blocks, is mandatory. Information from other files will only be used for
building blocks that come from an :ref:`rtp` file with the same base name. The
user can add building blocks to a force field by having additional files
with the same base name in their working directory. By default, only
extra building blocks can be defined, but calling :ref:`pdb2gmx <gmx pdb2gmx>` with the ``-rtpo``
option will allow building blocks in a local file to replace the default
ones in the force field.

Residue database
~~~~~~~~~~~~~~~~

The files holding the residue databases have the extension :ref:`rtp`.
Originally this file contained building blocks (amino acids) for
proteins, and is the |Gromacs| interpretation of the ``rt37c4.dat`` file of
GROMOS. So the residue database file contains information (bonds,
charges, charge groups, and improper dihedrals) for a frequently-used
building block. It is better *not* to change this file because it is
standard input for :ref:`pdb2gmx <gmx pdb2gmx>`, but if changes are needed make them in the
:ref:`top` file (see :ref:`topfile`), or in a :ref:`rtp` file in the working
directory as explained in sec. :ref:`pdb2gmxfiles`. Defining topologies
of new small molecules is probably easier by writing an include topology
file :ref:`itp` directly. This will be discussed in section :ref:`molitp`.
When adding a new protein residue to the database, don’t forget to add
the residue name to the residuetypes.dat file, so that :ref:`grompp <gmx grompp>`, :ref:`make_ndx <gmx make_ndx>`
and analysis tools can recognize the residue as a protein residue (see
:ref:`defaultgroups`).

The :ref:`rtp` files are only used by :ref:`pdb2gmx <gmx pdb2gmx>`. As mentioned before, the only
extra information this program needs from the :ref:`rtp` database is bonds,
charges of atoms, charge groups, and improper dihedrals, because the
rest is read from the coordinate input file. Some proteins contain
residues that are not standard, but are listed in the coordinate file.
You have to construct a building block for this “strange” residue,
otherwise you will not obtain a :ref:`top` file. This also holds for molecules
in the coordinate file such as ligands, polyatomic ions, crystallization
co-solvents, etc. The residue database is constructed in the following
way:

::

    [ bondedtypes ]  ; mandatory
    ; bonds  angles  dihedrals  impropers
         1       1          1          2  ; mandatory

    [ GLY ]  ; mandatory

     [ atoms ]  ; mandatory 
    ; name  type  charge  chargegroup 
         N     N  -0.280     0
         H     H   0.280     0
        CA   CH2   0.000     1
         C     C   0.380     2
         O     O  -0.380     2

     [ bonds ]  ; optional
    ;atom1 atom2      b0      kb
         N     H
         N    CA
        CA     C
         C     O
        -C     N

     [ exclusions ]  ; optional
    ;atom1 atom2

     [ angles ]  ; optional
    ;atom1 atom2 atom3    th0    cth

     [ dihedrals ]  ; optional
    ;atom1 atom2 atom3 atom4   phi0     cp   mult

     [ impropers ]  ; optional
    ;atom1 atom2 atom3 atom4     q0     cq
         N    -C    CA     H
        -C   -CA     N    -O

    [ ZN ]

     [ atoms ]
        ZN    ZN   2.000     0

The file is free format; the only restriction is that there can be at
most one entry on a line. The first field in the file is the ``[ bondedtypes ]`` field,
which is followed by four numbers, indicating the interaction type for
bonds, angles, dihedrals, and improper dihedrals. The file contains
residue entries, which consist of atoms and (optionally) bonds, angles,
dihedrals, and impropers. The charge group codes denote the charge group
numbers. Atoms in the same charge group should always be ordered
consecutively. When using the hydrogen database with :ref:`pdb2gmx <gmx pdb2gmx>` for adding
missing hydrogens (see :ref:`hdb`), the atom names defined in the :ref:`rtp`
entry should correspond exactly to the naming convention used in the
hydrogen database. The atom names in the bonded interaction can be
preceded by a minus or a plus, indicating that the atom is in the
preceding or following residue respectively. Explicit parameters added
to bonds, angles, dihedrals, and impropers override the standard
parameters in the :ref:`itp` files. This should only be used in special cases.
Instead of parameters, a string can be added for each bonded
interaction. This is used in GROMOS-96 :ref:`rtp` files. These strings are
copied to the topology file and can be replaced by force-field
parameters by the C-preprocessor in :ref:`grompp <gmx grompp>` using ``#define`` statements.

:ref:`pdb2gmx <gmx pdb2gmx>` automatically generates all angles. This means
that for most force fields the ``[ angles ]`` field is only
useful for overriding :ref:`itp` parameters. For the GROMOS-96
force field the interaction number of all angles needs to be specified.

:ref:`pdb2gmx <gmx pdb2gmx>` automatically generates one proper dihedral for every rotatable
bond, preferably on heavy atoms. When the ``[ dihedrals ]`` field is used, no other
dihedrals will be generated for the bonds corresponding to the specified
dihedrals. It is possible to put more than one dihedral function on a
rotatable bond. In the case of CHARMM27 FF :ref:`pdb2gmx <gmx pdb2gmx>` can add correction
maps to the dihedrals using the default ``-cmap`` option. Please refer to
:ref:`charmmff` for more information.

:ref:`pdb2gmx <gmx pdb2gmx>` sets the number of exclusions to 3, which means
that interactions between atoms connected by at most 3 bonds are
excluded. Pair interactions are generated for all pairs of atoms that
are separated by 3 bonds (except pairs of hydrogens). When more
interactions need to be excluded, or some pair interactions should not
be generated, an ``[ exclusions ]`` field can be added,
followed by pairs of atom names on separate lines. All non-bonded and
pair interactions between these atoms will be excluded.

Residue to building block database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each force field has its own naming convention for residues. Most
residues have consistent naming, but some, especially those with
different protonation states, can have many different names. The
:ref:`r2b` files are used to convert standard residue names to
the force-field build block names. If no :ref:`r2b` is present
in the force-field directory or a residue is not listed, the building
block name is assumed to be identical to the residue name. The
:ref:`r2b` can contain 2 or 5 columns. The 2-column format has
the residue name in the first column and the building block name in the
second. The 5-column format has 3 additional columns with the building
block for the residue occurring in the N-terminus, C-terminus and both
termini at the same time (single residue molecule). This is useful for,
for instance, the AMBER force fields. If one or more of the terminal
versions are not present, a dash should be entered in the corresponding
column.

There is a |Gromacs| naming convention for residues which is only apparent
(except for the :ref:`pdb2gmx <gmx pdb2gmx>` code) through the
:ref:`r2b` file and ``specbond.dat`` files. This
convention is only of importance when you are adding residue types to an
:ref:`rtp` file. The convention is listed in :numref:`Table %s <tab-r2b>`.
For special bonds with, for instance,
a heme group, the |Gromacs| naming convention is introduced through
``specbond.dat`` (see :ref:`specbond`),
which can subsequently be translated by the :ref:`r2b` file,
if required.

.. |NDEL| replace:: N\ :math:`_\delta`
.. |NEPS| replace:: N\ :math:`_\epsilon`

.. _tab-r2b:

.. table:: Internal |Gromacs| residue naming convention.

           +--------------+-----------------------------------------------------------+
           | |Gromacs| ID | Residue                                                   |
           +==============+===========================================================+
           | ARG          | protonated arginine                                       |
           +--------------+-----------------------------------------------------------+
           | ARGN         | neutral arginine                                          |
           +--------------+-----------------------------------------------------------+
           | ASP          | negatively charged aspartic acid                          |
           +--------------+-----------------------------------------------------------+
           | ASPH         | neutral aspartic acid                                     |
           +--------------+-----------------------------------------------------------+
           | CYS          | neutral cysteine                                          |
           +--------------+-----------------------------------------------------------+
           | CYS2         | cysteine with sulfur bound to another cysteine or a heme  |
           +--------------+-----------------------------------------------------------+
           | GLU          |  negatively charged glutamic acid                         |
           +--------------+-----------------------------------------------------------+
           | GLUH         |  neutral glutamic acid                                    |
           +--------------+------------------------------+----------------------------+
           | HISD         | neutral histidine with |NDEL| protonated                  |
           +--------------+-----------------------------------------------------------+
           | HISE         | neutral histidine with |NEPS| protonated                  |
           +--------------+------------------------------+----------------------------+
           | HISH         | positive histidine with both |NDEL| and |NEPS| protonated |
           +--------------+-----------------------------------------------------------+
           | HIS1         | histidine bound to a heme                                 |
           +--------------+-----------------------------------------------------------+
           | LYSN         | neutral lysine                                            |
           +--------------+-----------------------------------------------------------+
           | LYS          | protonated lysine                                         |
           +--------------+-----------------------------------------------------------+
           | HEME         | heme                                                      |
           +--------------+-----------------------------------------------------------+


Atom renaming database
~~~~~~~~~~~~~~~~~~~~~~

Force fields often use atom names that do not follow IUPAC or PDB
convention. The :ref:`arn` database is used to translate the
atom names in the coordinate file to the force-field names. Atoms that
are not listed keep their names. The file has three columns: the
building block name, the old atom name, and the new atom name,
respectively. The residue name supports question-mark wildcards that
match a single character.

An additional general atom renaming file called
``xlateat.dat`` is present in the ``share/top``
directory, which translates common non-standard atom names in the
coordinate file to IUPAC/PDB convention. Thus, when writing force-field
files, you can assume standard atom names and no further atom name
translation is required, except for translating from standard atom names
to the force-field ones.

Hydrogen database
~~~~~~~~~~~~~~~~~

The hydrogen database is stored in :ref:`hdb` files. It contains information
for the :ref:`pdb2gmx <gmx pdb2gmx>` program on how to connect hydrogen atoms to existing
atoms. In versions of the database before |Gromacs| 3.3, hydrogen atoms
were named after the atom they are connected to: the first letter of the
atom name was replaced by an ‘H.’ In the versions from 3.3 onwards, the
H atom has to be listed explicitly, because the old behavior was
protein-specific and hence could not be generalized to other molecules.
If more than one hydrogen atom is connected to the same atom, a number
will be added to the end of the hydrogen atom name. For example, adding
two hydrogen atoms to ``ND2`` (in asparagine), the hydrogen atoms will
be named ``HD21`` and ``HD22``. This is important since atom naming in
the :ref:`rtp` file (see :ref:`rtp`) must be the same. The format of the
hydrogen database is as follows:

::

    ; res   # additions
            # H add type    H       i       j       k
    ALA     1
            1       1       H       N       -C      CA
    ARG     4
            1       2       H       N       CA      C
            1       1       HE      NE      CD      CZ
            2       3       HH1     NH1     CZ      NE
            2       3       HH2     NH2     CZ      NE

On the first line we see the residue name (ALA or ARG) and the number of
kinds of hydrogen atoms that may be added to this residue by the
hydrogen database. After that follows one line for each addition, on
which we see:

-  The number of H atoms added

-  The method for adding H atoms, which can be any of:

   #. | *one planar hydrogen*, *e.g.* *rings or peptide bond*
      | One hydrogen atom (n) is generated, lying in the plane of atoms
        (i,j,k) on the plane bisecting angle (j-i-k) at a distance of
        0.1 nm from atom i, such that the angles (n-i-j) and (n-i-k) are
        :math:`>` 90\ :math:`^{\rm o}`.

   #. | *one single hydrogen*, *e.g.* *hydroxyl*
      | One hydrogen atom (n) is generated at a distance of 0.1 nm from
        atom i, such that angle (n-i-j)=109.5 degrees and dihedral
        (n-i-j-k)=trans.

   #. | *two planar hydrogens*, *e.g.* *ethylene -C=CH*:math:`_2`, *or amide
        -C(=O)NH*:math:`_2`
      | Two hydrogens (n1,n2) are generated at a distance of 0.1 nm from
        atom i, such that angle (n1-i-j)=(n2-i-j)=120 degrees and
        dihedral (n1-i-j-k)=cis and (n2-i-j-k)=trans, such that names
        are according to IUPAC standards \ :ref:`129 <refiupac70>`.

   #. | *two or three tetrahedral hydrogens*, *e.g.* *-CH*:math:`_3`
      | Three (n1,n2,n3) or two (n1,n2) hydrogens are generated at a
        distance of 0.1 nm from atom i, such that angle
        (n1-i-j)=(n2-i-j)=(n3-i-j)=109.47:math:`^{\rm o}`, dihedral
        (n1-i-j-k)=trans, (n2-i-j-k)=trans+120 and
        (n3-i-j-k)=trans+240:math:`^{\rm o}`.

   #. | *one tetrahedral hydrogen*, *e.g.* *C*\ :math:`_3`\* CH*
      | One hydrogen atom (n:math:`^\prime`) is generated at a distance
        of 0.1 nm from atom i in tetrahedral conformation such that
        angle
        (n:math:`^\prime`-i-j)=(n:math:`^\prime`-i-k)=(n:math:`^\prime`-i-l)=109.47:math:`^{\rm o}`.

   #. | *two tetrahedral hydrogens*, *e.g.* *C-CH*\ :math:`_2`\*-C*
      | Two hydrogen atoms (n1,n2) are generated at a distance of 0.1 nm
        from atom i in tetrahedral conformation on the plane bisecting
        angle j-i-k with angle
        (n1-i-n2)=(n1-i-j)=(n1-i-k)=109.47:math:`^{\rm o}`.

   #. | *two water hydrogens*
      | Two hydrogens are generated around atom i according to
        SPC \ :ref:`80 <refBerendsen81>` water geometry. The symmetry
        axis will alternate between three coordinate axes in both
        directions.

   #. | *three water “hydrogens”*
      | Two hydrogens are generated around atom i according to
        SPC \ :ref:`80 <refBerendsen81>` water geometry. The symmetry
        axis will alternate between three coordinate axes in both
        directions. In addition, an extra particle is generated on the
        position of the oxygen with the first letter of the name
        replaced by ‘M’. This is for use with four-atom water models
        such as TIP4P \ :ref:`128 <refJorgensen83>`.

   #. | *four water “hydrogens”*
      | Same as above, except that two additional particles are
        generated on the position of the oxygen, with names ‘LP1’ and
        ‘LP2.’ This is for use with five-atom water models such as
        TIP5P \ :ref:`130 <refMahoney2000a>`.

-  The name of the new H atom (or its prefix, *e.g.* ``HD2``
   for the asparagine example given earlier).

-  Three or four control atoms (i,j,k,l), where the first always is the
   atom to which the H atoms are connected. The other two or three
   depend on the code selected. For water, there is only one control
   atom.

Some more exotic cases can be approximately constructed from the above
tools, and with suitable use of energy minimization are good enough for
beginning MD simulations. For example secondary amine hydrogen, nitrenyl
hydrogen (:math:`\mathrm{C}=\mathrm{NH}`)
and even ethynyl hydrogen could be approximately constructed using
method 2 above for hydroxyl hydrogen.

Termini database
~~~~~~~~~~~~~~~~

The termini
databases
are stored in ``aminoacids.n.tdb`` and
``aminoacids.c.tdb`` for the N- and C-termini respectively.
They contain information for the :ref:`pdb2gmx <gmx pdb2gmx>` program on how
to connect new atoms to existing ones, which atoms should be removed or
changed, and which bonded interactions should be added. Their format is
as follows (from ``gromos43a1.ff/aminoacids.c.tdb``):

::

    [ None ]

    [ COO- ]
    [ replace ]
    C	C	C	12.011	0.27
    O 	O1	OM	15.9994	-0.635
    OXT	O2	OM	15.9994	-0.635
    [ add ]
    2	8	O	C	CA	N
    	OM	15.9994	-0.635
    [ bonds ]
    C	O1	gb_5
    C	O2	gb_5
    [ angles ]
    O1	C	O2	ga_37
    CA	C	O1	ga_21
    CA	C	O2	ga_21
    [ dihedrals ]
    N	CA	C	O2	gd_20
    [ impropers ]
    C	CA	O2	O1	gi_1

The file is organized in blocks, each with a header specifying the name
of the block. These blocks correspond to different types of termini that
can be added to a molecule. In this example ``[ COO- ]`` is
the first block, corresponding to changing the terminal carbon atom into
a deprotonated carboxyl group. ``[ None ]`` is the second
terminus type, corresponding to a terminus that leaves the molecule as
it is. Block names cannot be any of the following:
``replace``, ``add``, ``delete``,
``bonds``, ``angles``,
``dihedrals``, ``impropers``. Doing so would
interfere with the parameters of the block, and would probably also be
very confusing to human readers.

For each block the following options are present:

-  | ``[ replace ]``
   | Replace an existing atom by one with a different atom type, atom
     name, charge, and/or mass. This entry can be used to replace an
     atom that is present both in the input coordinates and in the
     :ref:`rtp` database, but also to only rename an atom in
     the input coordinates such that it matches the name in the force
     field. In the latter case, there should also be a corresponding
     ``[ add ]`` section present that gives instructions to
     add the same atom, such that the position in the sequence and the
     bonding is known. Such an atom can be present in the input
     coordinates and kept, or not present and constructed by
     :ref:`pdb2gmx <gmx pdb2gmx>`. For each atom to be replaced on line
     should be entered with the following fields:

   -  name of the atom to be replaced

   -  new atom name (optional)

   -  new atom type

   -  new mass

   -  new charge

-  | ``[ add ]``
   | Add new atoms. For each (group of) added atom(s), a two-line entry
     is necessary. The first line contains the same fields as an entry
     in the hydrogen database (name of the new atom, number of atoms,
     type of addition, control atoms, see :ref:`hdb`), but the
     possible types of addition are extended by two more, specifically
     for C-terminal additions:

   #. | *two carboxyl oxygens, -COO*:math:`^-`
      | Two oxygens (n1,n2) are generated according to rule 3, at a
        distance of 0.136 nm from atom i and an angle
        (n1-i-j)=(n2-i-j)=117 degrees

   #. | *carboxyl oxygens and hydrogen, -COOH*
      | Two oxygens (n1,n2) are generated according to rule 3, at
        distances of 0.123 nm and 0.125 nm from atom i for n1 and n2,
        respectively, and angles (n1-i-j)=121 and (n2-i-j)=115 degrees.
        One hydrogen (n:math:`^\prime`) is generated around n2 according
        to rule 2, where n-i-j and n-i-j-k should be read as
        n\ :math:`^\prime`-n2-i and n\ :math:`^\prime`-n2-i-j,
        respectively.

   After this line, another line follows that specifies the details of
   the added atom(s), in the same way as for replacing atoms, *i.e.*:

   -  atom type

   -  mass

   -  charge

   -  charge group (optional)

   Like in the hydrogen database (see :ref:`rtp`), when more than one
   atom is connected to an existing one, a number will be appended to
   the end of the atom name. **Note** that, like in the hydrogen
   database, the atom name is now on the same line as the control atoms,
   whereas it was at the beginning of the second line prior to |Gromacs|
   version 3.3. When the charge group field is left out, the added atom
   will have the same charge group number as the atom that it is bonded
   to.

-  | ``[ delete ]``
   | Delete existing atoms. One atom name per line.

-  | ``[ bonds ]``, ``[ angles ]``,
     ``[ dihedrals ]`` and ``[ impropers ]``
   | Add additional bonded parameters. The format is identical to that
     used in the :ref:`rtp` file, see :ref:`rtp`.

Virtual site database
~~~~~~~~~~~~~~~~~~~~~

Since we cannot rely on the positions of hydrogens in input files, we
need a special input file to decide the geometries and parameters with
which to add virtual site hydrogens. For more complex virtual site
constructs (*e.g.* when entire aromatic side chains are made rigid) we
also need information about the equilibrium bond lengths and angles for
all atoms in the side chain. This information is specified in the
:ref:`vsd` file for each force field. Just as for the termini,
there is one such file for each class of residues in the
:ref:`rtp` file.

The virtual site database is not really a very simple list of
information. The first couple of sections specify which mass centers
(typically called MCH\ :math:`_3`/MNH:math:`_3`) to use for
CH\ :math:`_3`, NH\ :math:`_3`, and NH\ :math:`_2` groups. Depending on
the equilibrium bond lengths and angles between the hydrogens and heavy
atoms we need to apply slightly different constraint distances between
these mass centers. **Note** that we do *not* have to specify the actual
parameters (that is automatic), just the type of mass center to use. To
accomplish this, there are three sections names ``[ CH3 ]``,
``[ NH3 ]``, and ``[ NH2 ]``. For each of these we expect three columns.
The first column is the atom type bound to the 2/3 hydrogens, the second
column is the next heavy atom type which this is bound, and the third
column the type of mass center to use. As a special case, in the
``[ NH2 ]`` section it is also possible to specify ``planar`` in the
second column, which will use a different construction without mass
center. There are currently different opinions in some force fields
whether an NH\ :math:`_2` group should be planar or not, but we try hard
to stick to the default equilibrium parameters of the force field.

The second part of the virtual site database contains explicit
equilibrium bond lengths and angles for pairs/triplets of atoms in
aromatic side chains. These entries are currently read by specific
routines in the virtual site generation code, so if you would like to
extend it *e.g.* to nucleic acids you would also need to write new code
there. These sections are named after the short amino acid names
(``[ PHE ]``, ``[ TYR ]``, ``[ TRP ]``, ``[ HID ]``, ``[ HIE ]``,
``[ HIP ]``), and simply contain 2 or 3 columns with atom names,
followed by a number specifying the bond length (in nm) or angle (in
degrees). **Note** that these are approximations of the equilibrated
geometry for the entire molecule, which might not be identical to the
equilibrium value for a single bond/angle if the molecule is strained.

.. _specbond:

Special bonds
~~~~~~~~~~~~~

The primary mechanism used by
:ref:`pdb2gmx <gmx pdb2gmx>` to generate
inter-residue bonds relies on head-to-tail linking of backbone atoms in
different residues to build a macromolecule. In some cases (*e.g.*
disulfide bonds, a heme
group, branched
polymers), it is necessary to
create inter-residue bonds that do not lie on the backbone. The file
``specbond.dat`` takes
care of this function. It is necessary that the residues belong to the
same ``[ moleculetype ]``. The ``-merge`` and
``-chainsep`` functions of :ref:`pdb2gmx <gmx pdb2gmx>` can be
useful when managing special inter-residue bonds between different
chains.

The first line of ``specbond.dat`` indicates the number of
entries that are in the file. If you add a new entry, be sure to
increment this number. The remaining lines in the file provide the
specifications for creating bonds. The format of the lines is as
follows:

``resA atomA nbondsA resB atomB nbondsB length newresA
newresB``

The columns indicate:

#. ``resA`` The name of residue A that participates in the
   bond.

#. ``atomA`` The name of the atom in residue A that forms
   the bond.

#. ``nbondsA`` The total number of bonds
   ``atomA`` can form.

#. ``resB`` The name of residue B that participates in the
   bond.

#. ``atomB`` The name of the atom in residue B that forms
   the bond.

#. ``nbondsB`` The total number of bonds
   ``atomB`` can form.

#. ``length`` The reference length for the bond. If
   ``atomA`` and ``atomB`` are not within
   ``length`` :math:`\pm` 10% in the coordinate file
   supplied to :ref:`pdb2gmx <gmx pdb2gmx>`, no bond will be formed.

#. ``newresA`` The new name of residue A, if necessary. Some
   force fields use *e.g.* CYS2 for a cysteine in a disulfide or heme
   linkage.

#. ``newresB`` The new name of residue B, likewise.

File formats
------------

.. _topfile:

Topology file
~~~~~~~~~~~~~

The topology file is built following the |Gromacs| specification for a
molecular topology. A :ref:`top` file can be generated by
:ref:`pdb2gmx <gmx pdb2gmx>`. All possible entries in the topology file are
listed in :numref:`Tables %s <tab-topfile1>` and
:numref:`%s <tab-topfile2>`. Also tabulated are: all the units of
the parameters, which interactions can be perturbed for free energy
calculations, which bonded interactions are used by
:ref:`grompp <gmx grompp>` for generating exclusions, and which bonded
interactions can be converted to constraints by :ref:`grompp <gmx grompp>`.

.. |VCR| replace:: V\ :math:`^{(cr)}`
.. |WCR| replace:: W\ :math:`^{(cr)}`
.. |CRO| replace:: :math:`^{(cr)}`
.. |TREF| replace:: :numref:`Table %s <tab-topfile2>`
.. |AKJM| replace:: :math:`a~\mathrm{kJ~mol}^{-1}`
.. |KJN6| replace:: :math:`\mathrm{kJ~mol}^{-1}~\mathrm{nm}^{-6}`
.. |BNM| replace:: :math:`b~\mathrm{nm}^{-1}`
.. |C6LJ| replace:: :math:`c_6`
.. |STAR| replace:: :math:`^{(*)}`
.. |NREX| replace:: :math:`n_{ex}^{(nrexcl)}`
.. |QEMU| replace:: :math:`q` (e); :math:`m` (u) 
.. |MQM| replace:: :math:`q,m`

.. _tab-topfile1:

.. table:: The :ref:`topology <top>` file.

        +------------------------------------------------------------------------------------------------------------+
        | Parameters                                                                                                 |
        +===================+===========================+=====+====+=========================================+=======+
        | interaction type  | directive                 | #   | f. | parameters                              | F. E. |
        |                   |                           | at. | tp |                                         |       |
        +-------------------+---------------------------+-----+----+-----------------------------------------+-------+
        | *mandatory*       | ``defaults``              |            non-bonded function type;                       |
        |                   |                           |            combination rule\ |CRO|;                        |
        |                   |                           |            generate pairs (no/yes);                        |
        |                   |                           |            fudge LJ (); fudge QQ ()                        |
        +-------------------+---------------------------+------------------------------------------------------------+
        | *mandatory*       | ``atomtypes``             |            atom type; m (u); q (e); particle type;         | 
        |                   |                           |            |VCR| ; |WCR|                                   |
        +-------------------+---------------------------+------------------------------------------------------------+
        |                   | ``bondtypes``             |  (see |TREF|, directive ``bonds``)                         |
        +                   +                           +                                                            +
        |                   | ``pairtypes``             |  (see |TREF|, directive ``pairs``)                         |
        +                   +                           +                                                            +
        |                   | ``angletypes``            |  (see |TREF|, directive ``angles``)                        |
        +                   +                           +                                                            +
        |                   | ``dihedraltypes``\ |STAR| |  (see |TREF|, directive ``dihedrals``)                     |
        +                   +                           +                                                            +
        |                   | ``constrainttypes``       |  (see |TREF|, directive ``constraints``)                   |
        +-------------------+---------------------------+-----+----+-------------------------------------------------+
        | LJ                | ``nonbond_params``        |  2  | 1  |  |VCR|  ; |WCR|                                 |
        +                   +                           +     +    +                                                 +
        | Buckingham        | ``nonbond_params``        |  2  | 2  |  |AKJM| ; |BNM|;                                |
        |                   |                           |     |    |  |C6LJ| (|KJN6|)                                |
        +-------------------+---------------------------+-----+----+-------------------------------------------------+

.. table:: 

        +------------------------------------------------------------------------------------------------------------+
        | Molecule definition(s)                                                                                     |
        +===================+===========================+============================================================+
        | *mandatory*       | ``moleculetype``          | molecule name; |NREX|                                      |
        +-------------------+---------------------------+-----+----------------------------------------------+-------+
        | *mandatory*       | ``atoms``                 | 1   | atom type; residue number;                   | type  |
        |                   |                           |     | residue name; atom name;                     |       |
        |                   |                           |     | charge group number; |QEMU|                  | |MQM| |
        +-------------------+---------------------------+-----+----------------------------------------------+-------+
        | intra-molecular interaction and geometry definitions as described in |TREF|                                |
        +------------------------------------------------------------------------------------------------------------+

.. table::

        +-------------+---------------+------------------------------------+
        | System      |               |                                    |
        +=============+===============+====================================+
        | *mandatory* | ``system``    | system name                        |
        +-------------+---------------+------------------------------------+
        | *mandatory* | ``molecules`` | molecule name; number of molecules |
        +-------------+---------------+------------------------------------+

.. table::

        +------------------------------+----------------------------------------------------+
        | Inter-molecular interactions |                                                    |
        +==============================+====================================================+
        | *optional*                   | ``intermolecular_interactions``                    |
        +------------------------------+----------------------------------------------------+
        | one or more bonded interactions as described in |TREF|, with two or more atoms,   |
        | no interactions that generate exclusions, no constraints, use global atom numbers |
        +-----------------------------------------------------------------------------------+

.. parsed-literal::

    '\# at' is the required number of atom type indices for this directive
    'f. tp' is the value used to select this function type
    'F. E.' indicates which of the parameters can be interpolated in free energy calculations
    |CRO| the combination rule determines the type of LJ parameters, see 
    |STAR| for ``dihedraltypes`` one can specify 4 atoms or the inner (outer for improper) 2 atoms
    |NREX| exclude neighbors :math:`n_{ex}` bonds away for non-bonded interactions
    For free energy calculations, type, :math:`q` and :math:`m`  or no parameters should be added
    for topology 'B' (:math:`\lambda = 1`) on the same line, after the normal parameters.

.. |BZERO| replace:: :math:`b_0`
.. |KB| replace:: :math:`k_b`
.. |KDR| replace:: :math:`k_{dr}`
.. |NM2| replace:: (kJ mol\ :math:`^{-1}`\ nm\ :math:`^{-2}`
.. |NM4| replace:: (kJ mol\ :math:`^{-1}`\ nm\ :math:`^{-4}`
.. |DKJ| replace:: :math:`D` (kJ mol\ :math:`^{-1}`
.. |BETA| replace:: :math:`\beta` (nm\ :math:`^{-1}`
.. |C23| replace:: :math:`C_{i=2,3}` (kJ mol\ :math:`^{-1}\ nm\ :math:`^{-i}`
.. |BMM| replace:: :math:`b_m`
.. |GE0| replace:: :math:`\geq 0`
.. |KO| replace:: :math:`k` 
.. |KJM| replace:: kJ mol\ :math:`^{-1}`
.. |LUU| replace:: low, up\ :math:`_1`,\ :math:`_2`
.. |MV| replace:: :math:`V`
.. |MW| replace:: :math:`W`
.. |QIJ| replace:: :math:`q_i`; :math:`q_j`
.. |THE0| replace:: :math:`\theta_0`
.. |KTHE| replace:: :math:`k_\theta`
.. |KJR2| replace:: kJ mol\ :math:`^{-1}`\ rad\ :math:`^{-2}`
.. |RN13| replace:: :math:`r_{13}`
.. |KUB| replace:: :math:`k_{UB}`
.. |C024| replace:: :math:`C_{i=0,1,2,3,4}`
.. |KJRI| replace:: kJ mol\ :math:`^{-1}`\ rad\ :math:`^{-i}`
.. |PHIS| replace:: :math:`\phi_s`
.. |PHI0| replace:: :math:`\phi_0`
.. |KPHI| replace:: :math:`k_\phi`
.. |PHIK| replace:: :math:`\phi,k`
.. |XI0| replace:: :math:`\xi_0`
.. |KXI| replace:: :math:`k_\xi`
.. |C0| replace:: :math:`C_0`
.. |C1| replace:: :math:`C_1`
.. |C2| replace:: :math:`C_2`
.. |C3| replace:: :math:`C_3`
.. |C4| replace:: :math:`C_4`
.. |C5| replace:: :math:`C_5`
.. |A0| replace:: :math:`a_0`
.. |A1| replace:: :math:`a_1`
.. |A2| replace:: :math:`a_2`
.. |A3| replace:: :math:`a_3`
.. |A4| replace:: :math:`a_4`
.. |DOH| replace:: :math:d_{\mbox{\sc oh}}`
.. |DHH| replace:: :math:d_{\mbox{\sc hh}}`
.. |AO| replace:: :math:`a`
.. |BO| replace:: :math:`b`
.. |CO| replace:: :math:`c`
.. |DO| replace:: :math:`d`
.. |KX| replace:: :math:`k_{x}`
.. |KY| replace:: :math:`k_{y}`
.. |KZ| replace:: :math:`k_{z}`
.. |GO| replace:: :math:`g`
.. |RO| replace:: :math:`r`
.. |DPHI| replace:: :math:`\Delta\phi`
.. |DIHR| replace:: :math:`k_{\mathrm{dihr}}`
.. |THET| replace:: :math:`\theta`
.. |NM| replace:: nm\ :math:`^{-1}`
.. |KC| replace:: :math:`k_c`
.. |THEK| replace:: :math:`\theta,k`
.. |R1E| replace:: :math:`r_{1e}`
.. |R2E| replace:: :math:`r_{2e}`
.. |R3E| replace:: :math:`r_{3e}`
.. |KRR| replace:: :math:`k_{rr'}`
.. |KRTH| replace:: :math:`k_{r\theta}`
.. |ALPH| replace:: :math:`\alpha`; |CO| (U nm\ :math:`^{\alpha}`
.. |UM1| replace:: U\ :math:`^{-1}`

.. _tab-topfile2:

.. table:: Details of ``[ moleculetype ]`` directives

            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | Name of interaction                | Topology file directive    | num.       | func.     | Order of parameters and their units                                     | use in     | 
            |                                    |                            | atoms [1]_ | type [2]_ |                                                                         | F.E.? [3]_ |
            +====================================+============================+============+===========+=========================================================================+============+
            | bond                               | ``bonds`` [4]_, [5]_       | 2          | 1         | |BZERO| (nm); |KB| |NM2|                                                | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | G96 bond                           | ``bonds`` [4]_, [5]_       | 2          | 2         | |BZERO| (nm); |KB| |NM4|                                                | all        |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | Morse                              | ``bonds`` [4]_, [5]_       | 2          | 3         | |BZERO| (nm); |DKJ|; |BETA|                                             | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | cubic bond                         | ``bonds`` [4]_, [5]_       | 2          | 4         | |BZERO| (nm); |C23|                                                     |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | connection                         | ``bonds`` [4]_             | 2          | 5         |                                                                         |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | harmonic potential                 | ``bonds``                  | 2          | 6         | |BZERO| (nm); |KB| |NM2|                                                | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | FENE bond                          | ``bonds`` [4]_             | 2          | 7         | |BMM|   (nm); |KB| |NM2|                                                |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | tabulated bond                     | ``bonds`` [4]_             | 2          | 8         | table number (|GE0|); |KO| |KJM|                                        | |KO|       |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | tabulated bond [6]_                | ``bonds``                  | 2          | 9         | table number (|GE0|); |KO| |KJM|                                        | |KO|       |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | restraint potential                | ``bonds``                  | 2          | 10        | |LUU| (nm); |KDR| (|NM2|)                                               | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | extra LJ or Coulomb                | ``pairs``                  | 2          | 1         | |MV| [7]_; |MW| [7]_                                                    | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | extra LJ or Coulomb                | ``pairs``                  | 2          | 2         | fudge QQ (); |QIJ| (e), |MV| [7]_; |MW| [7]_                            |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | extra LJ or Coulomb                | ``pairs_nb``               | 2          | 1         | |QIJ| (e); |MV| [7]_; |MW| [7]_                                         |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | angle                              | ``angles`` [5]_            | 3          | 1         | |THE0| (deg); |KTHE| (|KJR2|)                                           | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | G96 angle                          | ``angles`` [5]_            | 3          | 2         | |THE0| (deg); |KTHE| (|KJM|)                                            | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | cross bond-bond                    | ``angles``                 | 3          | 3         | |R1E|, |R2E| (nm); |KRR| (|NM2|)                                        |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | cross bond-angle                   | ``angles``                 | 3          | 4         | |R1E|, |R2E|, |R3E| (nm); |KRTH| (|NM2|)                                |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | Urey-Bradley                       | ``angles`` [5]_            | 3          | 5         | |THE0| (deg); |KTHE| (|KJR2|); |RN13| (nm); |KUB| (|NM2|)               | all        |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | quartic angle                      | ``angles`` [5]_            | 3          | 6         | |THE0| (deg); |C024| (|KJRI|)                                           |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | tabulated angle                    | ``angles``                 | 3          | 8         | table number (|GE0|); |KO| (|KJM|)                                      | |KO|       | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            |  |  restricted                     |                            |            |           |                                                                         |            |
            |  |  bending potential              | ``angles``                 | 3          | 10        | |THE0| (deg); |KTHE| (|KJM|)                                            |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | proper dihedral                    | ``dihedrals``              | 4          | 1         | |PHIS| (deg); |KPHI| (|KJM|); multiplicity                              | |PHIK|     | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | improper dihedral                  | ``dihedrals``              | 4          | 2         | |XI0| (deg); |KXI| (|KJR2|)                                             | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | Ryckaert-Bellemans dihedral        | ``dihedrals``              | 4          | 3         | |C0|, |C1|, |C2|, |C3|, |C4|, |C5| (|KJM|)                              | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | periodic improper dihedral         | ``dihedrals``              | 4          | 4         | |PHIS| (deg); |KPHI| (|KJM|); multiplicity                              | |PHIK|     | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | Fourier dihedral                   | ``dihedrals``              | 4          | 5         | |C1|, |C2|, |C3|, |C4|, |C5| (|KJM|)                                    | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | tabulated dihedral                 | ``dihedrals``              | 4          | 8         | table number (|GE0|); |KO| (|KJM|)                                      | |KO|       |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | proper dihedral (multiple)         | ``dihedrals``              | 4          | 9         | |PHIS| (deg); |KPHI| (|KJM|); multiplicity                              | |PHIK|     | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | restricted dihedral                | ``dihedrals``              | 4          | 10        | |PHI0| (deg); |KPHI| (|KJM|)                                            |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | combined bending-torsion potential | ``dihedrals``              | 4          | 11        | |A0|, |A1|, |A2|, |A3|, |A4| (|KJM|)                                    |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | exclusions                         | ``exclusions``             | 1          |           | one or more atom indices                                                |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | constraint                         | ``constraints`` [4]_       | 2          | 1         | |BZERO| (nm)                                                            | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | constraint [6]_                    | ``constraints``            | 2          | 2         | |BZERO| (nm)                                                            | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | SETTLE                             | ``settles``                | 1          | 1         | |DOH|, |DHH| (nm)                                                       |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | 2-body virtual site                | ``virtual_sites2``         | 3          | 1         | |AO| ()                                                                 |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | 3-body virtual site                | ``virtual_sites3``         | 4          | 1         | |AO|, |BO| ()                                                           |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | 3-body virtual site (fd)           | ``virtual_sites3``         | 4          | 2         | |AO| (); |DO| (nm)                                                      |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | 3-body virtual site (fad)          | ``virtual_sites3``         | 4          | 3         | |THET| (deg); |DO| (nm)                                                 |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | 3-body virtual site (out)          | ``virtual_sites3``         | 4          | 4         | |AO|, |BO| (); |CO| (|NM|)                                              |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | 4-body virtual site (fdn)          | ``virtual_sites4``         | 5          | 2         | |AO|, |BO| (); |CO| (nm)                                                |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | N-body virtual site (COG)          | ``virtual_sitesn``         | 1          | 1         | one or more constructing atom indices                                   |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | N-body virtual site (COM)          | ``virtual_sitesn``         | 1          | 2         | one or more constructing atom indices                                   |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | N-body virtual site (COW)          | ``virtual_sitesn``         | 1          | 3         |  |  one or more pairs consisting of                                     |            |
            |                                    |                            |            |           |  |  constructing atom index and weight                                  |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | position restraint                 | ``position_restraints``    | 1          | 1         | |KX|, |KY|, |KZ| (|NM2|)                                                | all        |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | flat-bottomed position restraint   | ``position_restraints``    | 1          | 2         | |GO|, |RO| (nm), |KO| (|NM2|)                                           |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | distance restraint                 | ``distance_restraints``    | 2          | 1         | type; label; |LUU| (nm); weight ()                                      |            | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | dihedral restraint                 | ``dihedral_restraints``    | 4          | 1         | |PHI0| (deg); |DPHI| (deg); |DIHR| (|KJR2|)                             | all        | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | orientation restraint              | ``orientation_restraints`` | 2          | 1         | exp.; label; |ALPH|; obs. (U); weight (|UM1|)                           |            |
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | angle restraint                    | ``angle_restraints``       | 4          | 1         | |THE0| (deg); |KC| (|KJM|); multiplicity                                | |THEK|     | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+
            | angle restraint (z)                | ``angle_restraints_z``     | 2          | 1         | |THE0| (deg); |KC| (|KJM|); multiplicity                                | |THEK|     | 
            +------------------------------------+----------------------------+------------+-----------+-------------------------------------------------------------------------+------------+

.. [1]
   The required number of atom indices for this directive
   
.. [2]
   The index to use to select this function type
   
.. [3]
   Indicates which of the parameters can be interpolated in free energy calculations
   
.. [4]
   This interaction type will be used by :ref:`grompp <gmx grompp>` for generating exclusions
   
.. [5]
   This interaction type can be converted to constraints by :ref:`grompp <gmx grompp>`
   
.. [7]
   The combination rule determines the type of LJ parameters, see
   
.. [6]
   No connection, and so no exclusions, are generated for this interaction

Description of the file layout:

-  Semicolon (;) and newline characters surround comments

-  On a line ending with :math:`\backslash` the newline character is
   ignored.

-  Directives are surrounded by ``[`` and ``]``

-  The topology hierarchy (which must be followed) consists of three
   levels:

   -  the parameter level, which defines certain force-field
      specifications (see :numref:`Table %s <tab-topfile1>`)

   -  the molecule level, which should contain one or more molecule
      definitions (see :numref:`Table %s <tab-topfile2>`)

   -  the system level, containing only system-specific information
      (``[ system ]`` and ``[ molecules ]``)

-  Items should be separated by spaces or tabs, not commas

-  Atoms in molecules should be numbered consecutively starting at 1

-  Atoms in the same charge group must be listed consecutively

-  The file is parsed only once, which implies that no forward
   references can be treated: items must be defined before they can be
   used

-  Exclusions can be generated from the bonds or overridden manually

-  The bonded force types can be generated from the atom types or
   overridden per bond

-  It is possible to apply multiple bonded interactions of the same type
   on the same atoms

-  Descriptive comment lines and empty lines are highly recommended

-  Starting with |Gromacs| version 3.1.3, all directives at the parameter
   level can be used multiple times and there are no restrictions on the
   order, except that an atom type needs to be defined before it can be
   used in other parameter definitions

-  If parameters for a certain interaction are defined multiple times
   for the same combination of atom types the last definition is used;
   starting with |Gromacs| version 3.1.3 :ref:`grompp <gmx grompp>` generates
   a warning for parameter redefinitions with different values

-  Using one of the ``[ atoms ]``,
   ``[ bonds ]``, ``[ pairs ]``,
   ``[ angles ]``, etc. without having used
   ``[ moleculetype ]`` before is meaningless and generates
   a warning

-  Using ``[ molecules ]`` without having used
   ``[ system ]`` before is meaningless and generates a
   warning.

-  After ``[ system ]`` the only allowed directive is
   ``[ molecules ]``

-  Using an unknown string in ``[ ]`` causes all the data
   until the next directive to be ignored and generates a warning

Here is an example of a topology file, ``urea.top``:

::

    ;
    ;       Example topology file
    ;
    ; The force-field files to be included
    #include "amber99.ff/forcefield.itp"

    [ moleculetype ]
    ; name  nrexcl
    Urea         3

    [ atoms ]
       1  C  1  URE      C      1     0.880229  12.01000   ; amber C  type
       2  O  1  URE      O      2    -0.613359  16.00000   ; amber O  type
       3  N  1  URE     N1      3    -0.923545  14.01000   ; amber N  type
       4  H  1  URE    H11      4     0.395055   1.00800   ; amber H  type
       5  H  1  URE    H12      5     0.395055   1.00800   ; amber H  type
       6  N  1  URE     N2      6    -0.923545  14.01000   ; amber N  type
       7  H  1  URE    H21      7     0.395055   1.00800   ; amber H  type
       8  H  1  URE    H22      8     0.395055   1.00800   ; amber H  type

    [ bonds ]
        1	2
        1	3	
        1   6
        3	4
        3	5
        6	7
        6	8

    [ dihedrals ] 
    ;   ai    aj    ak    al funct  definition
         2     1     3     4   9     
         2     1     3     5   9     
         2     1     6     7   9     
         2     1     6     8   9     
         3     1     6     7   9     
         3     1     6     8   9     
         6     1     3     4   9     
         6     1     3     5   9     

    [ dihedrals ] 
         3     6     1     2   4     
         1     4     3     5   4	 
         1     7     6     8   4

    [ position_restraints ]
    ; you wouldn't normally use this for a molecule like Urea,
    ; but we include it here for didactic purposes
    ; ai   funct    fc
       1     1     1000    1000    1000 ; Restrain to a point
       2     1     1000       0    1000 ; Restrain to a line (Y-axis)
       3     1     1000       0       0 ; Restrain to a plane (Y-Z-plane)

    [ dihedral_restraints ]
    ; ai   aj    ak    al  type  phi  dphi  fc
        3    6     1    2     1  180     0  10
        1    4     3    5     1  180     0  10

    ; Include TIP3P water topology
    #include "amber99/tip3p.itp"

    [ system ]
    Urea in Water

    [ molecules ]
    ;molecule name   nr.
    Urea             1
    SOL              1000

Here follows the explanatory text.

**#include “amber99.ff/forcefield.itp” :** this includes
the information for the force field you are using, including bonded and
non-bonded parameters. This example uses the AMBER99 force field, but
your simulation may use a different force field. :ref:`grompp <gmx grompp>`
will automatically go and find this file and copy-and-paste its content.
That content can be seen in
``share/top/amber99.ff/forcefield.itp}``, and it
is

::

    #define _FF_AMBER
    #define _FF_AMBER99

    [ defaults ]
    ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
    1               2               yes             0.5     0.8333

    #include "ffnonbonded.itp"
    #include "ffbonded.itp"

The two ``#define`` statements set up the conditions so that
future parts of the topology can know that the AMBER 99 force field is
in use.

**[ defaults ] :**

-  ``nbfunc`` is the non-bonded function type. Use 1 (Lennard-Jones) or 2
   (Buckingham)

-  ``comb-rule`` is the number of the combination rule (see :ref:`nbpar`).

-  ``gen-pairs`` is for pair generation. The default is
   ‘no’, *i.e.* get 1-4 parameters from the pairtypes list. When
   parameters are not present in the list, stop with a fatal error.
   Setting ‘yes’ generates 1-4 parameters that are not present in the
   pair list from normal Lennard-Jones parameters using
   ``fudgeLJ``

-  ``fudgeLJ`` is the factor by which to multiply
   Lennard-Jones 1-4 interactions, default 1

-  ``fudgeQQ`` is the factor by which to multiply
   electrostatic 1-4 interactions, default 1

-  :math:`N` is the power for the repulsion term in a 6-\ :math:`N`
   potential (with nonbonded-type Lennard-Jones only), starting with
   |Gromacs| version 4.5, :ref:`grompp <gmx mdrun>` also reads and applies
   :math:`N`, for values not equal to 12 tabulated interaction functions
   are used (in older version you would have to use user tabulated
   interactions).

**Note** that ``gen-pairs``, ``fudgeLJ``,
``fudgeQQ``, and :math:`N` are optional.
``fudgeLJ`` is only used when generate pairs is set to
‘yes’, and ``fudgeQQ`` is always used. However, if you want
to specify :math:`N` you need to give a value for the other parameters
as well.

Then some other ``#include`` statements add in the large
amount of data needed to describe the rest of the force field. We will
skip these and return to ``urea.top``. There we will see

**[ moleculetype ] :** defines the name of your molecule
in this :ref:`top` and nrexcl = 3 stands for excluding
non-bonded interactions between atoms that are no further than 3 bonds
away.

**[ atoms ] :** defines the molecule, where
``nr`` and ``type`` are fixed, the rest is user
defined. So ``atom`` can be named as you like,
``cgnr`` made larger or smaller (if possible, the total
charge of a charge group should be zero), and charges can be changed
here too.

**[ bonds ] :** no comment.

**[ pairs ] :** LJ and Coulomb 1-4 interactions

**[ angles ] :** no comment

**[ dihedrals ] :** in this case there are 9 proper
dihedrals (funct = 1), 3 improper (funct = 4) and no Ryckaert-Bellemans
type dihedrals. If you want to include Ryckaert-Bellemans type dihedrals
in a topology, do the following (in case of *e.g.* decane):

::

    [ dihedrals ]
    ;  ai    aj    ak    al funct       c0       c1       c2
        1    2     3     4     3 
        2    3     4     5     3

In the original implementation of the potential for
alkanes \ :ref:`131 <refRyckaert78>` no 1-4 interactions were used, which means that in
order to implement that particular force field you need to remove the
1-4 interactions from the ``[ pairs ]`` section of your
topology. In most modern force fields, like OPLS/AA or Amber the rules
are different, and the Ryckaert-Bellemans potential is used as a cosine
series in combination with 1-4 interactions.

**[ position_restraints ] :** harmonically restrain the selected particles to reference
positions (:ref:`positionrestraint`). The reference positions are read
from a separate coordinate file by :ref:`grompp <gmx grompp>`.

**[ dihedral_restraints ] :** restrain selected dihedrals to a reference value. The
implementation of dihedral restraints is described in section
:ref:`dihedralrestraint` of the manual. The parameters specified in
the ``[dihedral_restraints]`` directive are as follows:

-  ``type`` has only one possible value which is 1

-  ``phi`` is the value of :math:`\phi_0` in :eq:`eqn. %s <eqndphi>` and
   :eq:`eqn. %s <eqndihre>` of the manual.

-  ``dphi`` is the value of :math:`\Delta\phi` in :eq:`eqn. %s <eqndihre>` of the
   manual.

-  ``fc`` is the force constant :math:`k_{dihr}` in :eq:`eqn. %s <eqndihre>` of the
   manual.

**#include “tip3p.itp” :** includes a topology file that was already
constructed (see section :ref:`molitp`).

**[ system ] :** title of your system, user-defined

**[ molecules ] :** this defines the total number of (sub)molecules in your system
that are defined in this :ref:`top`. In this example file, it stands for 1
urea molecule dissolved in 1000 water molecules. The molecule type ``SOL``
is defined in the ``tip3p.itp`` file. Each name here must correspond to a
name given with ``[ moleculetype ]`` earlier in the topology. The order of the blocks of
molecule types and the numbers of such molecules must match the
coordinate file that accompanies the topology when supplied to :ref:`grompp <gmx grompp>`.
The blocks of molecules do not need to be contiguous, but some tools
(e.g. :ref:`genion <gmx genion>`) may act only on the first or last such block of a
particular molecule type. Also, these blocks have nothing to do with the
definition of groups (see sec. :ref:`groupconcept` and
sec. :ref:`usinggroups`).

.. _molitp:

Molecule.itp file
~~~~~~~~~~~~~~~~~

If you construct a topology file you will use frequently (like the water
molecule, ``tip3p.itp``, which is already constructed for
you) it is good to make a ``molecule.itp`` file. This only
lists the information of one particular molecule and allows you to
re-use the ``[ moleculetype ]`` in multiple systems without
re-invoking :ref:`pdb2gmx <gmx pdb2gmx>` or manually copying and pasting. An
example ``urea.itp`` follows:

::

    [ moleculetype ]
    ; molname	nrexcl
    URE		3

    [ atoms ]
       1  C  1  URE      C      1     0.880229  12.01000   ; amber C  type
    ...
       8  H  1  URE    H22      8     0.395055   1.00800   ; amber H  type

    [ bonds ]
        1	2
    ...
        6	8
    [ dihedrals ] 
    ;   ai    aj    ak    al funct  definition
         2     1     3     4   9     
    ...
         6     1     3     5   9     
    [ dihedrals ] 
         3     6     1     2   4     
         1     4     3     5   4	 
         1     7     6     8   4

Using :ref:`itp` files results in a very short
:ref:`top` file:

::

    ;
    ;       Example topology file
    ;
    ; The force field files to be included
    #include "amber99.ff/forcefield.itp"

    #include "urea.itp"

    ; Include TIP3P water topology
    #include "amber99/tip3p.itp"

    [ system ]
    Urea in Water

    [ molecules ]
    ;molecule name   nr.
    Urea             1
    SOL              1000

Ifdef statements
~~~~~~~~~~~~~~~~

A very powerful feature in |Gromacs| is the use of ``#ifdef``
statements in your :ref:`top` file. By making use of this
statement, and associated ``#define`` statements like were
seen in ``amber99.ff/forcefield.itp`` earlier,
different parameters for one molecule can be used in the same
:ref:`top` file. An example is given for TFE, where there is
an option to use different charges on the atoms: charges derived by De
Loof et al. :ref:`132 <refLoof92>` or by Van Buuren and
Berendsen \ :ref:`133 <refBuuren93a>`. In fact, you can use much of the
functionality of the C preprocessor, ``cpp``, because
:ref:`grompp <gmx grompp>` contains similar pre-processing functions to scan
the file. The way to make use of the ``#ifdef`` option is as
follows:

-  either use the option ``define = -DDeLoof`` in the
   :ref:`mdp` file (containing :ref:`grompp <gmx grompp>` input
   parameters), or use the line ``#define DeLoof`` early in
   your :ref:`top` or :ref:`itp` file; and

-  put the ``#ifdef`` statements in your
   :ref:`top`, as shown below:


::

    ...



    [ atoms ]
    ; nr     type     resnr    residu     atom      cgnr      charge        mass
    #ifdef DeLoof
    ; Use Charges from DeLoof
       1        C        1        TFE        C         1        0.74        
       2        F        1        TFE        F         1       -0.25        
       3        F        1        TFE        F         1       -0.25        
       4        F        1        TFE        F         1       -0.25        
       5      CH2        1        TFE      CH2         1        0.25        
       6       OA        1        TFE       OA         1       -0.65        
       7       HO        1        TFE       HO         1        0.41        
    #else
    ; Use Charges from VanBuuren
       1        C        1        TFE        C         1        0.59        
       2        F        1        TFE        F         1       -0.2         
       3        F        1        TFE        F         1       -0.2         
       4        F        1        TFE        F         1       -0.2         
       5      CH2        1        TFE      CH2         1        0.26        
       6       OA        1        TFE       OA         1       -0.55        
       7       HO        1        TFE       HO         1        0.3         
    #endif

    [ bonds ]
    ;  ai    aj funct           c0           c1
        6     7     1 1.000000e-01 3.138000e+05 
        1     2     1 1.360000e-01 4.184000e+05 
        1     3     1 1.360000e-01 4.184000e+05 
        1     4     1 1.360000e-01 4.184000e+05 
        1     5     1 1.530000e-01 3.347000e+05 
        5     6     1 1.430000e-01 3.347000e+05 
    ...

This mechanism is used by :ref:`pdb2gmx <gmx pdb2gmx>` to implement optional position
restraints (:ref:`positionrestraint`) by ``#include``-ing an :ref:`itp` file
whose contents will be meaningful only if a particular ``#define`` is set
(and spelled correctly!)

Topologies for free energy calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Free energy differences between two systems, A and B, can be calculated
as described in sec. :ref:`fecalc`. Systems A and B are described by
topologies consisting of the same number of molecules with the same
number of atoms. Masses and non-bonded interactions can be perturbed by
adding B parameters under the ``[ atoms ]`` directive. Bonded interactions can be
perturbed by adding B parameters to the bonded types or the bonded
interactions. The parameters that can be perturbed are listed in
:numref:`Tables %s <tab-topfile1>` and :numref:`%s <tab-topfile2>`.
The :math:`\lambda`-dependence of the
interactions is described in section sec. :ref:`feia`. The bonded
parameters that are used (on the line of the bonded interaction
definition, or the ones looked up on atom types in the bonded type
lists) is explained in :numref:`Table %s <tab-topfe>`. In most cases, things should
work intuitively. When the A and B atom types in a bonded interaction
are not all identical and parameters are not present for the B-state,
either on the line or in the bonded types, :ref:`grompp <gmx grompp>` uses the A-state
parameters and issues a warning. For free energy calculations, all or no
parameters for topology B (:math:`\lambda = 1`) should be added on the
same line, after the normal parameters, in the same order as the normal
parameters. From |Gromacs| 4.6 onward, if :math:`\lambda` is treated as a
vector, then the ``bonded-lambdas`` component controls all bonded terms that
are not explicitly labeled as restraints. Restrain terms are controlled
by the ``restraint-lambdas`` component.

.. |NOT| replace:: :math:`-`

.. _tab-topfe:

.. table:: The bonded parameters that are used for free energy topologies,
           on the line of the bonded interaction definition or looked up
           in the bond types section based on atom types. A and B indicate the
           parameters used for state A and B respectively, + and |NOT| indicate
           the (non-)presence of parameters in the topology, x indicates that
           the presence has no influence.

           +--------------------+---------------+---------------------------------+---------+
           | B-state atom types | parameters    | parameters in bonded types      |         |
           +                    +               +-----------------+---------------+         +
           | all identical to   | on line       | A atom types    | B atom types  | message |
           +                    +-------+-------+-------+---------+-------+-------+         +
           | A-state atom types | A     | B     | A     | B       | A     | B     |         |
           +====================+=======+=======+=======+=========+=======+=======+=========+
           |                    | +AB   | |NOT| | x     | x       |       |       |         |
           |                    | +A    | +B    | x     | x       |       |       |         |
           | yes                | |NOT| | |NOT| | |NOT| | |NOT|   |       |       | error   |
           |                    | |NOT| | |NOT| | +AB   | |NOT|   |       |       |         |
           |                    | |NOT| | |NOT| | +A    | +B      |       |       |         |
           +--------------------+-------+-------+-------+---------+-------+-------+---------+
           |                    | +AB   | |NOT| | x     | x       | x     | x     | warning |
           |                    | +A    | +B    | x     | x       | x     | x     |         |
           |                    | |NOT| | |NOT| | |NOT| | |NOT|   | x     | x     | error   |
           | no                 | |NOT| | |NOT| | +AB   | |NOT|   | |NOT| | |NOT| | warning |
           |                    | |NOT| | |NOT| | +A    | +B      | |NOT| | |NOT| | warning |
           |                    | |NOT| | |NOT| | +A    | x       | +B    | |NOT| |         |
           |                    | |NOT| | |NOT| | +A    | x       | +     | +B    |         |
           +--------------------+-------+-------+-------+---------+-------+-------+---------+



Below is an example of a topology which changes from 200 propanols to
200 pentanes using the GROMOS-96 force field.

::

     
    ; Include force field parameters
    #include "gromos43a1.ff/forcefield.itp"

    [ moleculetype ]
    ; Name            nrexcl
    PropPent          3

    [ atoms ]
    ; nr type resnr residue atom cgnr  charge    mass  typeB chargeB  massB
      1    H    1     PROP    PH    1   0.398    1.008  CH3     0.0  15.035
      2   OA    1     PROP    PO    1  -0.548  15.9994  CH2     0.0  14.027
      3  CH2    1     PROP   PC1    1   0.150   14.027  CH2     0.0  14.027
      4  CH2    1     PROP   PC2    2   0.000   14.027
      5  CH3    1     PROP   PC3    2   0.000   15.035

    [ bonds ]
    ;  ai    aj funct    par_A  par_B 
        1     2     2    gb_1   gb_26
        2     3     2    gb_17  gb_26
        3     4     2    gb_26  gb_26
        4     5     2    gb_26

    [ pairs ]
    ;  ai    aj funct
        1     4     1
        2     5     1

    [ angles ]
    ;  ai    aj    ak funct    par_A   par_B
        1     2     3     2    ga_11   ga_14
        2     3     4     2    ga_14   ga_14
        3     4     5     2    ga_14   ga_14

    [ dihedrals ]
    ;  ai    aj    ak    al funct    par_A   par_B
        1     2     3     4     1    gd_12   gd_17
        2     3     4     5     1    gd_17   gd_17

    [ system ]
    ; Name
    Propanol to Pentane

    [ molecules ]
    ; Compound        #mols
    PropPent          200

Atoms that are not perturbed, ``PC2`` and
``PC3``, do not need B-state parameter specifications, since
the B parameters will be copied from the A parameters. Bonded
interactions between atoms that are not perturbed do not need B
parameter specifications, as is the case for the last bond in the
example topology. Topologies using the OPLS/AA force field need no
bonded parameters at all, since both the A and B parameters are
determined by the atom types. Non-bonded interactions involving one or
two perturbed atoms use the free-energy perturbation functional forms.
Non-bonded interactions between two non-perturbed atoms use the normal
functional forms. This means that when, for instance, only the charge of
a particle is perturbed, its Lennard-Jones interactions will also be
affected when lambda is not equal to zero or one.

**Note** that this topology uses the GROMOS-96 force field, in which the
bonded interactions are not determined by the atom types. The bonded
interaction strings are converted by the C-preprocessor. The force-field
parameter files contain lines like:

::

    #define gb_26       0.1530  7.1500e+06

    #define gd_17     0.000       5.86          3

.. _constraintforce:

Constraint forces
~~~~~~~~~~~~~~~~~

| The constraint force between two atoms in one molecule can be
  calculated with the free energy perturbation code by adding a
  constraint between the two atoms, with a different length in the A and
  B topology. When the B length is 1 nm longer than the A length and
  lambda is kept constant at zero, the derivative of the Hamiltonian
  with respect to lambda is the constraint force. For constraints
  between molecules, the pull code can be used, see sec. :ref:`pull`.
  Below is an example for calculating the constraint force at 0.7 nm
  between two methanes in water, by combining the two methanes into one
  “molecule.” **Note** that the definition of a “molecule” in |Gromacs|
  does not necessarily correspond to the chemical definition of a
  molecule. In |Gromacs|, a “molecule” can be defined as any group of
  atoms that one wishes to consider simultaneously. The added constraint
  is of function type 2, which means that it is not used for generating
  exclusions (see sec. :ref:`excl`). Note that the constraint free energy
  term is included in the derivative term, and is specifically included
  in the ``bonded-lambdas`` component. However, the free energy for changing
  constraints is *not* included in the potential energy differences used
  for BAR and MBAR, as this requires reevaluating the energy at each of
  the constraint components. This functionality is planned for later
  versions.

::

    ; Include force-field parameters
    #include "gromos43a1.ff/forcefield.itp"

    [ moleculetype ]
    ; Name            nrexcl
    Methanes               1

    [ atoms ]
    ; nr   type   resnr  residu   atom    cgnr     charge    mass
       1    CH4     1     CH4      C1       1          0    16.043
       2    CH4     1     CH4      C2       2          0    16.043
    [ constraints ]
    ;  ai    aj funct   length_A  length_B
        1     2     2        0.7       1.7

    #include "gromos43a1.ff/spc.itp"

    [ system ]
    ; Name
    Methanes in Water

    [ molecules ]
    ; Compound        #mols
    Methanes              1
    SOL                2002

Coordinate file
~~~~~~~~~~~~~~~

Files with the :ref:`gro` file extension contain a molecular
structure in GROMOS-87 format. A sample piece is included below:

::

    MD of 2 waters, reformat step, PA aug-91
        6
        1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
        1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
        1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
        2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
        2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
        2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
       1.82060   1.82060   1.82060

This format is fixed, *i.e.* all columns are in a fixed position. If you
want to read such a file in your own program without using the |Gromacs|
libraries you can use the following formats:

**C-format:**
``“%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f”``

Or to be more precise, with title *etc.* it looks like this:

::

      "%s\n", Title
      "%5d\n", natoms
      for (i=0; (i<natoms); i++) {
        "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
          residuenr,residuename,atomname,atomnr,x,y,z,vx,vy,vz
      }
      "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
        box[X][X],box[Y][Y],box[Z][Z],
        box[X][Y],box[X][Z],box[Y][X],box[Y][Z],box[Z][X],box[Z][Y]

**Fortran format:**
``(i5,2a5,i5,3f8.3,3f8.4)``

So ``confin.gro`` is the |Gromacs| coordinate file and is
almost the same as the GROMOS-87 file (for GROMOS users: when used with
``ntx=7``). The only difference is the box for which |Gromacs|
uses a tensor, not a vector.

.. _fforganization:

Force field organization
------------------------

.. _fffiles:

Force-field files
~~~~~~~~~~~~~~~~~

Many force fields are available by default. Force fields are detected by
the presence of ``<name>.ff`` directories in the
``$GMXLIB/share/gromacs/top`` sub-directory and/or the
working directory. The information regarding the location of the force
field files is printed by :ref:`pdb2gmx <gmx pdb2gmx>` so you can easily keep
track of which version of a force field is being called, in case you
have made modifications in one location or another. The force fields
included with |Gromacs| are:

-  AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24,
   1999-2012, 2003)

-  AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)

-  AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29,
   461-469, 1996)

-  AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21,
   1049-1074, 2000)

-  AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65,
   712-725, 2006)

-  AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al.,
   Proteins 78, 1950-58, 2010)

-  AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)

-  CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)

-  GROMOS96 43a1 force field

-  GROMOS96 43a2 force field (improved alkane dihedrals)

-  GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)

-  GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)

-  GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)

-  GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856,
   DOI: 10.1007/s00249-011-0700-9)

-  OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

A force field is included at the beginning of a topology file with an
``#include`` statement followed by
``<name>.ff/forcefield.itp``. This statement includes the
force-field file, which, in turn, may include other force-field files.
All the force fields are organized in the same way. An example of the
``amber99.ff/forcefield.itp`` was shown in
:ref:`topfile`.

For each force field, there several files which are only used by
:ref:`pdb2gmx <gmx pdb2gmx>`. These are: residue databases
(:ref:`rtp`) the hydrogen
database (:ref:`hdb`), two
termini databases (``.n.tdb`` and ``.c.tdb``,
see ) and the atom type database
(:ref:`atp`), which
contains only the masses. Other optional files are described in sec. :ref:`pdb2gmxfiles`.

Changing force-field parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If one wants to change the parameters of few bonded interactions in a
molecule, this is most easily accomplished by typing the parameters
behind the definition of the bonded interaction directly in the
:ref:`top` file under the ``[ moleculetype ]``
section (see :ref:`topfile` for the format and units).
If one wants to change the parameters for all instances of a
certain interaction one can change them in the force-field file or add a
new ``[ ???types ]`` section after including the force
field. When parameters for a certain interaction are defined multiple
times, the last definition is used. As of |Gromacs| version 3.1.3, a
warning is generated when parameters are redefined with a different
value. Changing the Lennard-Jones parameters of an atom type is not
recommended, because in the GROMOS force fields the Lennard-Jones
parameters for several combinations of atom types are not generated
according to the standard combination rules. Such combinations (and
possibly others that do follow the combination rules) are defined in the
``[ nonbond_params ]`` section, and changing the
Lennard-Jones parameters of an atom type has no effect on these
combinations.

Adding atom types
~~~~~~~~~~~~~~~~~

As of |Gromacs| version 3.1.3, atom types can be added in an extra
``[ atomtypes ]`` section after the inclusion of the
normal force field. After the definition of the new atom type(s),
additional non-bonded and pair parameters can be defined. In pre-3.1.3
versions of |Gromacs|, the new atom types needed to be added in the
``[ atomtypes ]`` section of the force-field files, because
all non-bonded parameters above the last ``[ atomtypes ]``
section would be overwritten using the standard combination rules.
