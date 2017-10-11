Common errors when using |Gromacs|
==================================

This page covers some of the errors that might be encountered when using |Gromacs|. The list
is by no means complete, but should help with investigating a number of issues.

.. Again, text shamelessly stolen from the old webpage

.. _common-errors:

Common errors during usage
--------------------------

.. _out-of-memory:

Out of memory when allocating
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The executed script has attempted to assign memory to be used in the calculation, but is unable to due to insufficient memory.

Possible solutions are:

* install more memory in the computer.
* use a computer with more memory.
* reduce the scope of the number of atoms selected for analysis.
* reduce the length of trajectory file being processed.
* in some cases confusion between Ångström and nm may lead to users wanting to generate a 
  :ref:`pdb2gmx <gmx pdb2gmx>` water box that is |10to3| times larger than what they think it is (e.g. genbox).

.. |10to3| replace:: 10\ :sup:`3`

The user should bear in mind that the cost in time and/or memory for various activities will 
scale with the number of atoms/groups/residues *N* or the simulation length *T* as order N,
NlogN, or |Nsquared| (or maybe worse!) and the same for *T*, depending on the type of activity.
If it takes a long time, have a think about what you are doing, and the underlying algorithm
(see the `Reference manual`_, man page, or use the -h flag for the utility), and
see if there's something sensible you can do that has better scaling properties.

.. _Reference manual: `gmx-manual-parent-dir`_ 
.. |Nsquared| replace:: N\ :sup:`2`

.. _pdb2gmx-errors:

Errors in pdb2gmx
-----------------

Residue 'XXX' not found in residue topology database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


This means that the force field you have selected while running :ref:`pdb2gmx <gmx pdb2gmx>` does not have an entry in
the :ref:`residue database<rtp>` for XXX. The :ref:`residue database<rtp>` entry is necessary both for stand-alone
molecules (e.g. formaldehyde) or a peptide (standard or non-standard). This entry defines the atom
types, connectivity, bonded and non-bonded interaction types for the residue and is necessary
to use :ref:`pdb2gmx <gmx pdb2gmx>` to build a :ref:`top` file. A :ref:`residue database<rtp>`
entry may be missing simply because the
database does not contain the residue at all, or because the name is different.

For new users, this error appears because they are running :ref:`pdb2gmx <gmx pdb2gmx>` blindly on a
:ref:`PDB<pdb>` file they have without consideration of the contents of the file. A :ref:`force field<gmx-force-field>`
is not something that is magical, it can only deal with molecules or residues (building blocks) that are
provided in the :ref:`residue database<rtp>` or included otherwise.

If you want to use :ref:`pdb2gmx <gmx pdb2gmx>` to automatically generate your topology, you have
to ensure that the appropriate :ref:`rtp` entry is present within the desired :ref:`force field<gmx-force-field>` and
has the same name as the building block you are trying to use. If you call your
molecule "HIS," then :ref:`pdb2gmx <gmx pdb2gmx>` will not magically build a random molecule; it will try to
build histidine, based on the [ HIS ] entry in the :ref:`rtp` file, so it will look for the exact atomic entries for histidine, no more no less.

If you want a :ref:`topology<top>` for an arbitrary molecule, you cannot use :ref:`pdb2gmx <gmx pdb2gmx>` (unless you
build the :ref:`rtp` entry yourself). You will have to build it by hand, or use another program
(such as :ref:`x2top<gmx x2top>` or one of the scripts contributed by users) to build the :ref:`top` file.

If there is not an entry for this :ref:`residue<gmx-residue>` in the database, then
the options for obtaining the force field parameters are:

* see if there is a different name being used for the :ref:`residue<gmx-residue>` in the :ref:`residue database<rtp>` and rename as appropriate,
* parameterize the residue / molecule yourself (lots of work, even for an expert),
* find a :ref:`topology file<top>` for the molecule, convert it to an :ref:`itp` file and include it in your :ref:`top` file,
* use another :ref:`force field<gmx-force-field>` which has parameters available for this,
* search the primary literature for publications for parameters for the residue that are consistent with the force field that is being used.

Once you have determined the parameters and topology for your residue, see adding a residue to a force field for instructions on how to proceed.

Long bonds and/or missing atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are probably atoms missing earlier in the :ref:`pdb` file which makes :ref:`pdb2gmx<gmx pdb2gmx>` go crazy.
Check the screen output of :ref:`pdb2gmx<gmx pdb2gmx>`, as it will tell you which one is missing. Then add
the atoms in your :ref:`pdb` file, :ref:`energy minimization<gmx-energy-min>` will put them in the right place, or
fix the side chain with e.g. the :ref:`WhatIF` program.

Chain identifier 'X' was used in two non-sequential blocks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WARNING: atom X is missing in residue XXX Y in the pdb file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Atom X in residue YYY not found in rtp entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No force fields found (files with name 'forcefield.itp' in subdirectories ending on '.ff')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
