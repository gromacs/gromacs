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

.. raw:: latex

    \clearpage


