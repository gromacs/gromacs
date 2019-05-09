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

   #. | *two planar hydrogens*, *e.g.* *ethylene -C=CH*\ :math:`_2`, *or amide
        -C(=O)NH*\ :math:`_2`
      | Two hydrogens (n1,n2) are generated at a distance of 0.1 nm from
        atom i, such that angle (n1-i-j)=(n2-i-j)=120 degrees and
        dihedral (n1-i-j-k)=cis and (n2-i-j-k)=trans, such that names
        are according to IUPAC standards \ :ref:`129 <refiupac70>`.

   #. | *two or three tetrahedral hydrogens*, *e.g.* *-CH*\ :math:`_3`
      | Three (n1,n2,n3) or two (n1,n2) hydrogens are generated at a
        distance of 0.1 nm from atom i, such that angle
        (n1-i-j)=(n2-i-j)=(n3-i-j)=109.47\ :math:`^{\rm o}`, dihedral
        (n1-i-j-k)=trans, (n2-i-j-k)=trans+120 and
        (n3-i-j-k)=trans+240\ :math:`^{\rm o}`.

   #. | *one tetrahedral hydrogen*, *e.g.* *C*\ :math:`_3`\ *CH*
      | One hydrogen atom (n\ :math:`^\prime`) is generated at a distance
        of 0.1 nm from atom i in tetrahedral conformation such that
        angle
        (n\ :math:`^\prime`-i-j)=(n\ :math:`^\prime`-i-k)=(n\ :math:`^\prime`-i-l)=109.47\ :math:`^{\rm o}`.

   #. | *two tetrahedral hydrogens*, *e.g.* *C-CH*\ :math:`_2`\ *-C*
      | Two hydrogen atoms (n1,n2) are generated at a distance of 0.1 nm
        from atom i in tetrahedral conformation on the plane bisecting
        angle j-i-k with angle
        (n1-i-n2)=(n1-i-j)=(n1-i-k)=109.47\ :math:`^{\rm o}`.

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
        One hydrogen (n\ :math:`^\prime`) is generated around n2 according
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
(typically called MCH\ :math:`_3`/MNH\ :math:`_3`) to use for
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
