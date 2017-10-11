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

Errors in :ref:`pdb2gmx <gmx pdb2gmx>`
--------------------------------------

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

* see if there is a different name being used for the :ref:`residue <gmx-residue>` in the :ref:`residue database<rtp>` and rename as appropriate,
* parameterize the residue / molecule yourself (lots of work, even for an expert),
* find a :ref:`topology file<top>` for the molecule, convert it to an :ref:`itp` file and include it in your :ref:`top` file,
* use another :ref:`force field<gmx-force-field>` which has parameters available for this,
* search the primary literature for publications for parameters for the residue that are consistent with the force field that is being used.

Once you have determined the parameters and topology for your residue, see adding a residue to a force field for instructions on how to proceed.

Long bonds and/or missing atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are probably atoms missing earlier in the :ref:`pdb` file which makes :ref:`pdb2gmx <gmx pdb2gmx>` go crazy.
Check the screen output of :ref:`pdb2gmx <gmx pdb2gmx>`, as it will tell you which one is missing. Then add
the atoms in your :ref:`pdb` file, :ref:`energy minimization <gmx-energy-min>` will put them in the right place, or
fix the side chain with e.g. the `WHAT IF <http://swift.cmbi.ru.nl/whatif/>`_ program.


Chain identifier 'X' was used in two non-sequential blocks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This means that within the :ref:`coordinate file<gmx-structure-files>` fed to :ref:`pdb2gmx<gmx pdb2gmx>`, the X
chain has been split, possibly by the incorrect insertion of one molecule within another.
The solution is simple: move the inserted molecule to a location within the file so that it is not splitting another molecule.
This message may also mean that the same chain identifier has been used for two 
separate chains. In that case, rename the second chain to a unique identifier.

WARNING: atom X is missing in residue XXX Y in the pdb file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Related to the long bonds/missing atoms error above, this error is usually quite 
obvious in its meaning. That is, :ref:`pdb2gmx<gmx pdb2gmx>` expects certain atoms within 
the given residue, based on the entries in the force field :ref:`rtp` file. There are several cases to which this error applies:

* Missing hydrogen atoms; the error message may be suggesting that an entry in the :ref:`hdb`
  file is missing.  More likely, the nomenclature of your hydrogen atoms simply does not match
  what is expected by the :ref:`rtp` entry.  In this case, use ``-ignh`` to
  allow :ref:`pdb2gmx<gmx pdb2gmx>` to add the correct hydrogens for you,
  or re-name the problematic atoms.
* A terminal residue (usually the N-terminus) is missing H atoms; this usually suggests 
  that the proper ``-ter`` option has not been supplied or chosen properly. In the case of
  the :ref:`AMBER force fields<gmx-amber-ff>`, nomenclature is typically the problem.
  N-terminal and C-terminal residues must be prefixed by N and C, respectively.
  For example, an N-terminal alanine should not be listed in the :ref:`pdb` file
  as ALA, but rather NALA, as specified in the `ffamber <http://ffamber.cnsm.csulb.edu/ffamber.php>`_ instructions.
* Atoms are simply missing in the structure file provided to :ref:`pdb2gmx<gmx pdb2gmx>`;
  look for REMARK 465 and REMARK 470 entries in the :ref:`pdb` file. These atoms
  will have to be modeled in using external software. There is no GROMACS tool to re-construct incomplete models.

Contrary to what the error message says, the use of the option ``-missing``
is almost always inappropriate.  The ``-missing`` option should only be used to
generate specialized topologies for amino acid-like molecules to take 
advantage of :ref:`rtp` entries.  If you find yourself using ``-missing``
in order to generate a topology for a protein or nucleic acid,
don't; the topology produced is likely physically unrealistic.

Atom X in residue YYY not found in rtp entry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are attempting to assemble a topology using :ref:`pdb2gmx <gmx pdb2gmx>`, the atom names
are expected to match those found in the :ref:`rtp` file that define the building
block(s) in your structure.  In most cases, the problem arises from a naming mismatch,
so simply re-name the atoms in your :ref:`coordinate file <gmx-structure-files>` appropriately.
In other cases, you may be supplying a structure that has residues that do not conform
to the expectations of the `force field <gmx-force-field>`, in which case you should
investigate why such a difference is occurring and make a decision based on what you
find - use a different `force field <gmx-force-field>`, manually edit the structure, etc.

No force fields found (files with name 'forcefield.itp' in subdirectories ending on '.ff')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This means your environment is not configured to use GROMACS properly, because
:ref:`pdb2gmx <gmx pdb2gmx>` cannot find its databases of forcefield information. This could
happen because a GROMACS installation was moved from one location to another.
Either follow the instructions about
:ref:`getting access to GROMACS after installation <getting access to |Gromacs|>`
or re-install GROMACS before doing so.

Errors in :ref:`grompp <gmx grompp>`
------------------------------------

Found a second defaults directive file 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is caused by the ``[defaults]`` directive appearing more than once in the :ref:`topology <top>` or
:ref:`force field <gmx-force-field>` files for the system - it can only appear once. A typical cause of
this is a second defaults being set in an included :ref:`topology <top>` file, :ref:`itp`, that
has been sourced from somewhere else. For specifications on how the topology files work,
see the `reference manual`_, Section 5.6.::

    [ defaults ]
    ; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
    1       1       no       1.0       1.0

One solution is to simply comment out (or delete) the lines of code out in the file where it is included for the second time i.e.,::

    ;[ defaults ]
    ; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
    ;1       1       no       1.0       1.0

A better approach to finding a solution is to re-think what you are doing. The ``[defaults]``
directive should only be appearing at the top of your :ref:`top` file
where you choose the :ref:`force field <gmx-force-field>`. If you are trying
to mix two :ref:`force fields <gmx-force-field>`, then you are asking for trouble.
If a molecule :ref:`itp` file tries to choose a force field, then whoever produced it is asking for trouble.

Invalid order for directive xxx
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The directives in the .top and .itp files have rules about the order in which they can
appear, and this error is seen when the order is violated. Consider the examples and
discussion in chapter 5 of the `reference manual`_, and/or from tutorial material.
The :ref:`include file mechanism <gmx-topo-include>` cannot be used to ``#include`` a
file in just any old location, because they contain directives and these have to be properly placed.

In particular, ``Invalid order for directive defaults`` is a result of defaults being
set in the :ref:`topology <top>` or :ref:`force field <gmx-force-field>` files in the inappropriate location;
the ``[defaults]`` section can only appear once and must be the first directive in
the :ref:`topology <top>`.  The ``[defaults]`` directive is typically present in the :ref:`force field <gmx-force-field>`
file (forcefield.itp), and is added to the :ref:`topology <top>` when you ``#include`` this file in the system topology.

If the directive in question is ``[atomtypes]`` (which is the most common source of this error) or
any other bonded or nonbonded ``[*types]`` directive, typically the user is adding some
non-standard species (ligand, solvent, etc) that introduces new atom types or parameters
into the system. As indicated above, these new types and parameters must appear before
any ``[moleculetype]`` directive. The :ref:`force field <gmx-force-field>` has to be
fully constructed before any molecules can be defined.

Atom index n in position_restraints out of bounds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A common problem is placing position restraint files for multiple molecules out of order.
Recall that a position restraint :ref:`itp` file containing a ``[ position_restraints ]``
block can only belong to the ``[ moleculetype ]`` block that contains it. For example:

WRONG::

    #include "topol_A.itp"
    #include "topol_B.itp"
    #include "ligand.itp"
    #ifdef POSRES
    #include "posre_A.itp"
    #include "posre_B.itp"
    #include "ligand_posre.itp"
    #endif

RIGHT::

    #include "topol_A.itp"
    #ifdef POSRES
    #include "posre_A.itp"
    #endif
    #include "topol_B.itp"
    #ifdef POSRES
    #include "posre_B.itp"
    #endif
    #include "ligand.itp"
    #ifdef POSRES
    #include "ligand_posre.itp"
    #endif

Further, the atom index of each ``[position_restraint]`` must be relative to the
``[moleculetype]``, not relative to the system (because the parsing has not reached
``[molecules]`` yet, there is no such concept as "system"). So you cannot use the output 
of a tool like :ref:`genrestr <gmx genrestr>` blindly (as ``genrestr -h`` warns).

System has non-zero total charge 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Notifies you that counter-ions may be required for the system to neutralize the charge or
there may be problems with the topology.

If the charge is a non-integer, then this indicates that there is a problem with the :ref:`topology <top>`.
If :ref:`pdb2gmx <gmx pdb2gmx>` has been used, then look at the right hand comment column of the atom listing, which lists
the cumulative charge. This should be an integer after every residue (and/or charge group where
applicable). This will assist in finding the residue where things start departing from
integer values. Also check the capping groups that have been used.

If the charge is already close to an integer, then the difference is caused by
:ref:`rounding errors <gmx-floating-point>` and not a major problem.

Note for :ref:`PME <gmx-PME>` users: It is possible to use a uniform neutralizing background
charge in :ref:`PME <gmx-PME>` to compensate for a system with a net background charge. 
There is probably nothing wrong with this
in principle, because the uniform charge will not perturb the dynamics. Nevertheless, it is
standard practice to actually add counter-ions to make the system net neutral.

Incorrect number of parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Look at the :ref:`topology <top>` file for the system. You've not given enough parameters for one of the
bonded definitions.  Sometimes this also occurs if you've mangled the :ref:`Include File Mechanism <gmx-topo-include>`
or the topology file format (see: `reference manual`_ Chapter 5) when you edited the file.

Number of coordinates in coordinate file does not match topology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is pointing out that, based on the information provided in the :ref:`topology <top>` file, :ref:`top`,
the total number of atoms or particles within the system does not match exactly with what
is provided within the :ref:`coordinate file <gmx-structure-files>`, often a :ref:`gro` or a :ref:`pdb`.

The most common reason for this is simply that the user has failed to update the topology file
after solvating or adding additional molecules to the system, or made a typographical error in
the number of one of the molecules within the system. Ensure that the end of the topology file
being used contains something like the following, that matches exactly with what is within the
coordinate file being used, in terms of both numbers and order of the molecules::

    [ molecules ]
    ; Compound   #mol
    Protein      1
    SOL          10189
    NA+          10

Fatal error: No such moleculetype XXX
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each type of molecule in your ``[ molecules ]`` section of your :ref:`top` file must have a
corresponding ``[ moleculetype ]`` section defined previously, either in the :ref:`top` file or
an :ref:`included <gmx-topo-include>` :ref:`itp` file. See the `reference manual`_ section 5.6.1
for the syntax description. Your :ref:`top` file doesn't have such a definition for the
indicated molecule. Check the contents of the relevant files, how you have named your
molecules, and how you have tried to refer to them later. Pay attention to the status
of ``#ifdef`` and / or ``#include`` statements.

T-Coupling group XXX has fewer than 10% of the atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to specify separate :ref:`thermostats <gmx-thermostats>` (temperature coupling groups)
for every molecule type within a simulation. This is a particularly bad practice employed by
many new users to Molecular Dynamics Simulations.  Doing so is a bad idea, as you can
introduce errors and artifacts that are hard to predict. In some cases it is best to have all
molecules within a single group, using system. If separate coupling groups are required to avoid
the ``hot solvent cold solute`` problem, then ensure that they are of ``sufficient size`` and
combine molecule types that appear together within the simulation. For example, for
a protein in water with counter-ions, one would likely want to use ``Protein`` and ``Non-Protein``.

The cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element. Increase the box size or decrease rlist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error is generated in the cases as noted within the message. The dimensions of the box are such that an atom will
interact with itself (when using periodic boundary conditions), thus violating the minimum image convention.
Such an event is totally unrealistic and will introduce some serious artefacts. The solution is again what is
noted within the message, either increase the size of the simulation box so that it is at an absolute minimum
twice the cut-off length in all three dimensions (take care here if are using pressure coupling,
as the box dimensions will change over time and if they decrease even slightly, you will still be
violating the minimum image convention) or decrease the cut-off length (depending on the
:ref:`force field <gmx-force-field>` utilised, this may not be an option).

Atom index (1) in bonds out of bounds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This kind of error looks like::

    Fatal error:
    [ file spc.itp, line 32 ]
    Atom index (1) in bonds out of bounds (1-0).
    This probably means that you have inserted topology
    section "settles" in a part belonging to a different 
    molecule than you intended to. in that case move the
    "settles" section to the right molecule.

This error is fairly self-explanatory. You should look at your :ref:`top` file and check that all
of the ``[molecules]`` sections contain all of the data pertaining to that molecule, and no
other data. That is, you cannot ``#include`` another molecule type (:ref:`itp` file) before
the previous ``[moleculetype]`` has ended.Consult the examples in chapter 5 of the `reference manual`_
for information on the required ordering of the different ``[sections]``. Pay attention to
the contents of any files you have :ref:`included <gmx-topo-include>` with ``#include`` directives.

This error can also arise if you are using a water model that is not enabled for use with your
chosen :ref:`force field <gmx-force-field>` by default. For example, if you are attempting to use
the SPC water model with an :ref:`AMBER force field <gmx-amber-ff>`, you will see this error.
The reason is that, in ``spc.itp``, there is no ``#ifdef`` statement defining atom types for any
of the :ref:`AMBER force fields <gmx-amber-ff>`. You can either add this section yourself, or use a different water model.

XXX non-matching atom names
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error usually indicates that the order of the :ref:`topology <top>` file does not match that
of the :ref:`coordinate file <gmx-structure-files>`.  When running :ref:`grompp <gmx grompp>`, the
program reads through the :ref:`topology <top>`, mapping the supplied parameters to the atoms in
the :ref:`coordinate <gmx-structure-files>` file.  If there is a mismatch, this error is generated.
To remedy the problem, make sure that the contents of your ``[ molecules ]`` directive
matches the exact order of the atoms in the coordinate file.  

In some cases, the error is harmless. For example, when running simulations with the
`MARTINI force field <http://cgmartini.nl/>`_, the workflow relies on :ref:`grompp <gmx grompp>` to apply the
correct names, which are not previously assigned.  Also, perhaps you are using a
:ref:`coordinate <gmx-structure-files>` file that has the old (pre-4.5) ion nomenclature.
In this case, allowing :ref:`grompp <gmx grompp>` to re-assign names is harmless.
For just about any other situation, when this error comes up, **it should not be ignored**.
Just because the ``-maxwarn`` option is available does not mean you should use it in the blind
hope of your simulation working. It will undoubtedly :ref:`blow up <blowing-up>`.

The sum of the two largest charge group radii (X) is larger than rlist - rvdw/rcoulomb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error warns that some combination of settings will result in poor energy conservation at the
longest cutoff, which occurs when charge groups move in or out of neighborlist range.
The error can have two sources:

* Your charge groups encompass too many atoms. Most charge groups should be less than 4 atoms or less.
* Your :ref:`mdp` settings are incompatible with the chosen algorithms. For switch or shift functions,
  rlist must be larger than the longest cutoff (rvdw or rcoulomb) to provide buffer space for charge
  groups that move beyond the neighbor searching radius. If set incorrectly, you may miss
  interactions, contributing to poor energy conservation.

A similar error ("The sum of the two largest charge group radii (X) is larger than rlist") can arise under two circumstances:

* The charge groups are inappropriately large or rlist is set too low.
* Molecules are broken across periodic boundaries, which is not a problem in a periodic system.
  In this case, the sum of the two largest charge groups will correspond to a value of twice
  the box vector along which the molecule is broken.



