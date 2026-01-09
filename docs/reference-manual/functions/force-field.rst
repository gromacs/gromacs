Force field
-----------

A force field is built up from two distinct components:

-  The set of equations (called the *potential functions*) used to
   generate the potential energies and their derivatives, the forces.
   These are described in detail in the previous chapter.

-  The parameters used in this set of equations. These are not given in
   this manual, but in the data files corresponding to your |Gromacs|
   distribution.

Within one set of equations various sets of parameters can be used. Care
must be taken that the combination of equations and parameters form a
consistent set. It is in general dangerous to make *ad hoc* changes in a
subset of parameters, because the various contributions to the total
force are usually interdependent. This means in principle that every
change should be documented, verified by comparison to experimental data
and published in a peer-reviewed journal before it can be used.

|Gromacs| |version| includes several force fields, and
additional ones are available on the website. If you do not know which
one to select we recommend GROMOS-96 for united-atom setups and
OPLS-AA/L for all-atom parameters. That said, we describe the available
options in some detail.

GROMOS-96
~~~~~~~~~

.. warning::
    The GROMOS force fields have been parametrized with a physically
    incorrect multiple-time-stepping scheme for a twin-range cut-off. When
    used with a single-range cut-off (or a correct Trotter
    multiple-time-stepping scheme), physical properties, such as the density,
    might differ from the intended values. Since there are researchers
    actively working on validating GROMOS with modern integrators we have not
    yet removed the GROMOS force fields, but you should be aware of these
    issues and check if molecules in your system are affected before
    proceeding. Further information is available in
    `GitLab Issue 2884 <https://gitlab.com/gromacs/gromacs/-/issues/2884>`_ ,
    and a longer explanation of our
    decision to remove physically incorrect algorithms can be found at
    `DOI:10.26434/chemrxiv.11474583.v1 <https://doi.org/10.26434/chemrxiv.11474583.v1>`_ .

|Gromacs| supports the GROMOS-96 force fields \ :ref:`77 <refgromos96>`. All
parameters for the 43A1, 43A2 (development, improved alkane dihedrals),
45A3, 53A5, and 53A6 parameter sets are included. All standard building
blocks are included and topologies can be built automatically by
:ref:`pdb2gmx <gmx pdb2gmx>`.

The GROMOS-96 force field is a further development of the GROMOS-87
force field. It has improvements over the GROMOS-87 force field for
proteins and small molecules. **Note** that the sugar parameters present
in 53A6 do correspond to those published in 2004\ :ref:`110 <refOostenbrink2004>`,
which are different from those present in 45A4, which is not
included in |Gromacs| at this time. The 45A4 parameter set corresponds to
a later revision of these parameters. The GROMOS-96 force field is not,
however, recommended for use with long alkanes and lipids. The GROMOS-96
force field differs from the GROMOS-87 force field in a few respects:

-  the force field parameters

-  the parameters for the bonded interactions are not linked to atom
   types

-  a fourth power bond stretching potential (:ref:`G96bond`)

-  an angle potential based on the cosine of the angle
   (:ref:`G96angle`)

There are two differences in implementation between |Gromacs| and
GROMOS-96 which can lead to slightly different results when simulating
the same system with both packages:

-  in GROMOS-96 neighbor searching for solvents is performed on the
   first atom of the solvent molecule. This is not implemented in
   |Gromacs|, but the difference with searching by centers of charge
   groups is very small

-  the virial in GROMOS-96 is molecule-based. This is not implemented in
   |Gromacs|, which uses atomic virials

The GROMOS-96 force field was parameterized with a Lennard-Jones cut-off
of 1.4 nm, so be sure to use a Lennard-Jones cut-off
(``rvdw``) of at least 1.4. A larger cut-off is possible
because the Lennard-Jones potential and forces are almost zero beyond
1.4 nm.

GROMOS-96 files
^^^^^^^^^^^^^^^

|Gromacs| can read and write GROMOS-96 coordinate and trajectory files.
These files should have the extension :ref:`g96`. Such a file
can be a GROMOS-96 initial/final configuration file, a coordinate
trajectory file, or a combination of both. The file is fixed format; all
floats are written as 15.9, and as such, files can get huge. |Gromacs|
supports the following data blocks in the given order:

-  Header block:

   ::

       TITLE (mandatory)

-  Frame blocks:

   ::

       TIMESTEP (optional)
       POSITION/POSITIONRED (mandatory)
       VELOCITY/VELOCITYRED (optional)
       BOX (optional)

See the GROMOS-96 manual \ :ref:`77 <refgromos96>` for a complete
description of the blocks. **Note** that all |Gromacs| programs can read
compressed (.Z) or gzipped (.gz) files.

OPLS/AA
~~~~~~~

AMBER
~~~~~

|Gromacs| provides native support for the following AMBER force fields:

-  AMBER94 \ :ref:`111 <refCornell1995>`

-  AMBER96 \ :ref:`112 <refKollman1996>`

-  AMBER99 \ :ref:`113 <refWang2000>`

-  AMBER99SB \ :ref:`114 <refHornak2006>`

-  AMBER99SB-ILDN \ :ref:`115 <refLindorff2010>`

-  AMBER03 \ :ref:`116 <refDuan2003>`

-  AMBERGS \ :ref:`117 <refGarcia2002>`

-  AMBER14SB \ :ref:`197 <refMaier2015>`

-  AMBER19SB \ :ref:`198 <refTian2020>`

Dihedrals
^^^^^^^^^

|Gromacs| supports reordering atoms in select dihedrals to match
topologies produced by AMBER LEaP.

When this atom reordering is enabled, dihedral atoms
:math:`a_i`, :math:`a_j`, :math:`a_k`, and :math:`a_l` are reordered to
:math:`(a_l, a_k, a_j, a_i)` instead of :math:`(a_i, a_j, a_k, a_l)` if
:math:`at_i > at_l` or :math:`at_i = at_l` and :math:`at_j > at_k`, where
operator greater than is considered in lexicographical ordering,
:math:`at_i`, :math:`at_j`, :math:`at_k`, and :math:`at_l` are atom type
strings from the matching dihedral type, except in case of the wildcard
type X, which is replaced by space. This reordering doesn't affect the
resulting energies, but it is a necessity for what follows.

Additionally, but only for improper dihedrals, the third atom (:math:`a_k`
or :math:`a_j`) is taken to be central, while the remaining atoms
(:math:`a_i`, :math:`a_j`, and :math:`a_l` or :math:`a_l`, :math:`a_k`,
and :math:`a_i`) are then reordered by atom type from the matching
improper dihedral type. This can affect the resulting energies.

The reordering is performed by :ref:`gmx grompp` only when the
preprocessor define ``_FF_AMBER_LEAP_ATOM_REORDERING`` is encountered,
which enables this functionality only for select force fields (such as
AMBER19SB). In addition, preprocessor define-based enablement allows
fine-grained control in the topology file by combining ``#define`` and
``#undef`` statements, which can be useful for further force field
research and experimentation.

CMAPs
^^^^^

AMBER19SB and newer versions provide support for amino-acid-specific
energy correction maps (CMAPs). When these force fields are used with
:ref:`pdb2gmx <gmx pdb2gmx>`, the default option is to enable CMAPs.
If that is not the desired behavior, they can be disabled with ``-nocmap``
parameters.

CMAP types are specified using the following format:

::

    [ cmaptypes ]
    C-* N-GLY XC-GLY C-GLY N-* 1 24 24\
    3.44619344 4.59474328 4.73235504 5.80145072 8.84401368 ...

where:

- first five columns are pairs of atom and residue types,
- sixth coulumn is function type (presently only ``1`` is supported),
- seventh and eighth coulumns are grid size in ``x`` and ``y`` dimensions,
- following ``grid size x * grid size y`` values constitute the CMAP grid.

Note that AMBER19SB and newer versions require residue names to express
the residue type so that CMAP angle terms can be correctly assigned in
:ref:`grompp <gmx grompp>`. This is handled by the default of
``-rtpres auto`` for that option of :ref:`pdb2gmx <gmx pdb2gmx>` in
concert with the ``[ bondedtypes ]`` for the residue types for such
forcefields.

Note also that in the current implementation the grid size must be
identical for ``x`` and ``y`` dimensions and the same for all CMAP type
entries.

.. _charmmff:

CHARMM
~~~~~~

|Gromacs| supports the CHARMM force field for
proteins \ :ref:`118 <refmackerell04>`, :ref:`119 <refmackerell98>`,
lipids \ :ref:`120 <reffeller00>` and nucleic
acids \ :ref:`121 <reffoloppe00>`, :ref:`122 <refMac2000>`. The protein
parameters (and to some extent
the lipid and nucleic acid parameters) were thoroughly tested – both by
comparing potential energies between the port and the standard parameter
set in the CHARMM molecular simulation package, as well by how the
protein force field behaves together with |Gromacs|-specific techniques
such as virtual sites (enabling long time steps) recently
implemented \ :ref:`123 <refLarsson10>` – and the details and results are
presented in the paper by Bjelkmar et al. \ :ref:`124 <refBjelkmar10>`.
The nucleic acid parameters, as well as the ones for HEME, were
converted and tested by Michel Cuendet.

When selecting the CHARMM force field in
:ref:`pdb2gmx <gmx pdb2gmx>` the default option
is to use CMAP (for torsional correction map).
To exclude CMAP, use ``-nocmap``. The basic form of the CMAP
term implemented in |Gromacs| is a function of the :math:`\phi` and
:math:`\psi` backbone torsion angles. This term is defined in the
``rtp`` file by a ``[ cmap ]`` statement at the
end of each residue supporting CMAP. The following five atom names
define the two torsional angles. Atoms 1-4 define :math:`\phi`, and
atoms 2-5 define :math:`\psi`. The corresponding atom types are then
matched to the correct CMAP type in the ``cmap.itp`` file
that contains the correction maps.

A port of the CHARMM36 force field for use with |Gromacs| is also
available at `the MacKerell lab webpage <http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs>`_.

For branched polymers or other topologies not supported by
:ref:`pdb2gmx <gmx pdb2gmx>`, it is possible to
use TopoTools \ :ref:`125 <refkohlmeyer2016>` to generate a |Gromacs| top
file.

.. _cgforcefields:

Coarse-grained force fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coarse-graining is a systematic way of reducing the
number of degrees of freedom representing a system of interest. To
achieve this, typically whole groups of atoms are represented by single
beads and the coarse-grained force fields describes their effective
interactions. Depending on the choice of parameterization, the
functional form of such an interaction can be complicated and often
tabulated potentials are used.

Coarse-grained models are designed to reproduce certain properties of a
reference system. This can be either a full atomistic model or even
experimental data. Depending on the properties to reproduce there are
different methods to derive such force fields. An incomplete list of
methods is given below:

-  Conserving free energies

   -  Simplex method

   -  MARTINI force field (see next section)

-  Conserving distributions (like the radial distribution function),
   so-called structure-based coarse-graining

   -  (iterative) Boltzmann inversion

   -  Inverse Monte Carlo

-  Conversing forces

   -  Force matching

Note that coarse-grained potentials are state dependent (e.g.
temperature, density,...) and should be re-parametrized depending on the
system of interest and the simulation conditions. This can for example
be done using the Versatile Object-oriented Toolkit for Coarse-Graining
Applications (VOTCA) (**???**). The package was designed to assists in
systematic coarse-graining, provides implementations for most of the
algorithms mentioned above and has a well tested interface to |Gromacs|.
It is available as open source and further information can be found at
`www.votca.org <http://www.votca.org>`_.

MARTINI
~~~~~~~

The `MARTINI force field <https://cgmartini.nl>`__ is a coarse-grain parameter set
that allows for the construction of many systems, including proteins and membranes.

PLUM
~~~~

The PLUM force field :ref:`126 <refbereau12>` is an example of a solvent-free
protein-membrane model for which the membrane was derived from
structure-based coarse-graining \ :ref:`127 <refwang_jpcb10>`. A |Gromacs|
implementation can be found at
`github.com/tbereau/plumx <https://github.com/tbereau/plumx>`__.

