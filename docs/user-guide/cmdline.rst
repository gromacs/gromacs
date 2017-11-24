Command-line reference
======================

.. toctree::
   :hidden:
   :glob:

   /onlinehelp/gmx
   /onlinehelp/gmx-*

|Gromacs| includes many tools for preparing, running and analyzing
molecular dynamics simulations. These are all structured as part of a single
:command:`gmx` wrapper binary, and invoked with commands like :command:`gmx grompp`.
:ref:`mdrun <gmx mdrun>` is the only other binary that
:ref:`can be built <building just the mdrun binary>`; in the normal
build it can be run with :command:`gmx mdrun`. Documentation for these can
be found at the respective sections below, as well as on man pages (e.g.,
:manpage:`gmx-grompp(1)`) and with :samp:`gmx help {command}` or
:samp:`gmx {command} -h`.

If you've installed an MPI version of |Gromacs|, by default the
:command:`gmx` binary is called :command:`gmx_mpi` and you should adapt
accordingly.

Command-line interface and conventions
--------------------------------------

All |Gromacs| commands require an option before any arguments (i.e., all
command-line arguments need to be preceded by an argument starting with a
dash, and values not starting with a dash are arguments to the preceding
option).  Most options, except for boolean flags, expect an argument (or
multiple in some cases) after the option name.
The argument must be a separate command-line argument, i.e., separated by
space, as in ``-f traj.xtc``.  If more than one argument needs to be given to
an option, they should be similarly separated from each other.
Some options also have default arguments, i.e., just specifying the option
without any argument uses the default argument.
If an option is not specified at all, a default value is used; in the case of
optional files, the default might be not to use that file (see below).

All |Gromacs| command options start with a single dash, whether they are
single- or multiple-letter options.  However, two dashes are also recognized
(starting from 5.1).

In addition to command-specific options, some options are handled by the
:command:`gmx` wrapper, and can be specified for any command.  See
:doc:`wrapper binary help </onlinehelp/gmx>` for the list of such options.
These options are recognized both before the command name (e.g.,
:command:`gmx -quiet grompp`) as well as after the command name (e.g.,
:command:`gmx grompp -quiet`).
There is also a ``-hidden`` option that can be specified in combination with
``-h`` to show help for advanced/developer-targeted options.

Most analysis commands can process a trajectory with fewer atoms than the
run input or structure file, but only if the trajectory consists of the
first *n* atoms of the run input or structure file.

Handling specific types of command-line options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

boolean options
  Boolean flags can be specified like ``-pbc`` and negated like ``-nopbc``.
  It is also possible to use an explicit value like ``-pbc no`` and
  ``-pbc yes``.
file name options
  Options that accept files names have features that support using default file
  names (where the default file name is specific to that option):

  * If a required option is not set, the default is used.
  * If an option is marked optional, the file is not used unless the option
    is set (or other conditions make the file required).
  * If an option is set, and no file name is provided, the default is used.

  All such options will accept file names without a file extension.
  The extension is automatically appended in such a case.
  When multiple input formats are accepted, such as a generic structure format,
  the directory will be searched for files of each type with the supplied or
  default name. When no file with a recognized extension is found, an error is given.
  For output files with multiple formats, a default file type will be used.

  Some file formats can also be read from compressed (:file:`.Z` or
  :file:`.gz`) formats.
enum options
  Enumerated options (enum) should be used with one of the arguments listed in
  the option description. The argument may be abbreviated, and the first match
  to the shortest argument in the list will be selected.
vector options
  Some options accept a vector of values.  Either 1 or 3 parameters can be
  supplied; when only one parameter is supplied the two other values are also
  set to this value.
selection options
  See :doc:`/onlinehelp/selections`.

Commands by name
----------------

.. include:: /fragments/byname.rst

Commands by topic
-----------------

.. include:: /fragments/bytopic.rst

Special topics
--------------

The information in these topics is also accessible through
:samp:`gmx help {topic}` on the command line.

Selection syntax and usage
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::

   /onlinehelp/selections

.. _command-changes:

Command changes between versions
--------------------------------

Starting from |Gromacs| 5.0, some of the analysis commands (and a few other
commands as well) have changed significantly.

One main driver for this has been that many new tools mentioned below now
accept selections through one or more command-line options instead of prompting
for a static index group.  To take full advantage of selections, the interface
to the commands has changed somewhat, and some previous command-line options
are no longer present as the same effect can be achieved with suitable
selections.
Please see :doc:`/onlinehelp/selections` additional information on how to use
selections.

In the process, some old analysis commands have been removed in favor of more
powerful functionality that is available through an alternative tool.
For removed or replaced commands, this page documents how to perform the same
tasks with new tools.
For new commands, a brief note on the available features is given.  See the
linked help for the new commands for a full description.

This section lists only major changes; minor changes like additional/removed
options or bug fixes are not typically included.

Version 2018
^^^^^^^^^^^^

gmx trajectory
..............

**new**

:ref:`gmx trajectory` has been introduced as a selection-enabled version of
:ref:`gmx traj`.  It supports output of coordinates, velocities, and/or forces
for positions calculated for selections.

Version 2016
^^^^^^^^^^^^

Analysis on arbitrary subsets of atoms
......................................

Tools implemented in the new analysis framework can now operate upon trajectories
that match only a subset of the atoms in the input structure file.

gmx insert-molecules
....................

**improved**

:ref:`gmx insert-molecules` has gained an option ``-replace`` that makes it
possible to insert molecules into a solvated configuration, replacing any
overlapping solvent atoms.  In a fully solvated box, it is also possible to
insert into a certain region of the solvent only by selecting a subset of the
solvent atoms (``-replace`` takes a selection that can also contain expressions
like ``not within 1 of ...``).

gmx rdf
.......

**improved**

The normalization for the output RDF can now also be the radial number density.

gmx genconf
...........

**simplified**

Removed ``-block``, ``-sort`` and ``-shuffle``.

Version 5.1
^^^^^^^^^^^

General
.......

Symbolic links from 5.0 are no longer supported.  The only way to invoke a
command is through :samp:`gmx {<command>}`.

gmx pairdist
............

**new**

:ref:`gmx pairdist` has been introduced as a selection-enabled replacement for
:ref:`gmx mindist` (``gmx mindist`` still exists unchanged).  It can calculate
min/max pairwise distances between a pair of selections, including, e.g.,
per-residue minimum distances or distances from a single point to a set of
residue-centers-of-mass.

gmx rdf
.......

**rewritten**

:ref:`gmx rdf` has been rewritten for 5.1 to use selections for specifying the
points from which the RDFs are calculated.  The interface is mostly the same,
except that there are new command-line options to specify the selections.
The following additional changes have been made:

* ``-com`` and ``-rdf`` options have been removed.  Equivalent functionality is
  available through selections:

  * ``-com`` can be replaced with a :samp:`com of {<selection>}` as the
    reference selection.
  * ``-rdf`` can be replaced with a suitable set of selections (e.g.,
    :samp:`res_com of {<selection>}`) and/or using ``-seltype``.

* ``-rmax`` option is added to specify a cutoff for the RDFs.  If set to a
  value that is significantly smaller than half the box size, it can speed up
  the calculation significantly if a grid-based neighborhood search can be
  used.
* ``-hq`` and ``-fade`` options have been removed, as they are simply
  postprocessing steps on the raw numbers that can be easily done after the
  analysis.

Version 5.0
^^^^^^^^^^^

General
.......

Version 5.0 introduced the :command:`gmx` wrapper binary.
For backwards compatibility, this version still creates symbolic links by default for
old tools: e.g., ``g_order <options>`` is equivalent to ``gmx order <options>``, and
``g_order`` is simply a symbolic link on the file system.

g_bond
......

**replaced**

This tool has been removed in 5.0. A replacement is :ref:`gmx distance`.

You can provide your existing index file to :ref:`gmx distance`, and it will
calculate the same distances.  The differences are:

* ``-blen`` and ``-tol`` options have different default values.
* You can control the output histogram with ``-binw``.
* ``-aver`` and ``-averdist`` options are not present.  Instead, you can choose
  between the different things to calculate using ``-oav`` (corresponds to
  ``-d`` with ``-averdist``), ``-oall`` (corresponds to ``-d`` without
  ``-averdist``), ``-oh`` (corresponds to ``-o`` with ``-aver``), and
  ``-oallstat`` (corresponds to ``-l`` without ``-aver``).

You can produce any combination of output files.  Compared to ``g_bond``,
``gmx distance -oall`` is currently missing labels for the output columns.

g_dist
......

**replaced**

This tool has been removed in 5.0.  A replacement is :ref:`gmx distance` (for
most options) or :ref:`gmx select` (for ``-dist`` or ``-lt``).

If you had index groups A and B in :file:`index.ndx` for ``g_dist``, you can use the
following command to compute the same distance with ``gmx distance``::

  gmx distance -n index.ndx -select 'com of group "A" plus com of group "B"' -oxyz -oall

The ``-intra`` switch is replaced with ``-nopbc``.

If you used ``-dist D``, you can do the same calculation with ``gmx select``::

  gmx select -n index.ndx -select 'group "B" and within D of com of group "A"' -on/-oi/-os/-olt

You can select the output option that best suits your post-processing needs
(``-olt`` is a replacement for ``g_dist -dist -lt``)

gmx distance
............

**new**

:ref:`gmx distance` has been introduced as a selection-enabled replacement for
various tools that computed distances between fixed pairs of atoms (or
centers-of-mass of groups).  It has a combination of the features of ``g_bond``
and ``g_dist``, allowing computation of one or multiple distances, either
between atom-atom pairs or centers-of-mass of groups, and providing a
combination of output options that were available in one of the tools.

gmx gangle
..........

**new**

:ref:`gmx gangle` has been introduced as a selection-enabled replacement for
``g_sgangle``.  In addition to supporting atom-atom vectors, centers-of-mass
can be used as endpoints of the vectors, and there are a few additional angle
types that can be calculated.  The command also has basic support for
calculating normal angles between three atoms and/or centers-of-mass, making it
a partial replacement for :ref:`gmx angle` as well.

gmx protonate
.............

**replaced**

This was a very old tool originally written for united atom force fields,
where it was necessary to generate all hydrogens after running a trajectory
in order to calculate e.g. distance restraint violations. The functionality
to simply protonate a structure is available in :ref:`gmx pdb2gmx`. 
If there is significant interest, we might reintroduce it after moving to new
topology formats in the future.

gmx freevolume
..............

**new**

This tool has been introduced in 5.0.  It uses a Monte Carlo sampling method to
calculate the fraction of free volume within the box (using a probe of a given
size).

g_sas
.....

**rewritten**

This tool has been rewritten in 5.0, and renamed to :ref:`gmx sasa` (the
underlying surface area calculation algorithm is still the same).

The main difference in the new tool is support for selections.  Instead of
prompting for an index group, a (potentially dynamic) selection for the
calculation can be given with ``-surface``.  Any number of output groups can be
given with ``-output``, allowing multiple parts of the surface area to be
computed in a single run.  The total area of the ``-surface`` group is now
always calculated.

The tool no longer automatically divides the surface into hydrophobic and
hydrophilic areas, and there is no ``-f_index`` option.  The same effects can
be obtained by defining suitable selections for ``-output``.  If you want
output that contains the same numbers as with the old tool for a calculation
group ``A`` and output group ``B``, you can use ::

  gmx sasa -surface 'group "A"' -output '"Hydrophobic" group "A" and charge {-0.2 to 0.2}; "Hydrophilic" group "B" and not charge {-0.2 to 0.2}; "Total" group "B"'

Solvation free energy estimates are now calculated only if separately requested
with ``-odg``, and are written into a separate file.

Output option ``-i`` for a position restraint file is not currently implemented
in the new tool, but would not be very difficult to add if requested.

g_sgangle
.........

**replaced**

This tool has been removed in 5.0.  A replacement is :ref:`gmx gangle` (for
angle calculation) and :ref:`gmx distance` (for ``-od``, ``-od1``, ``-od2``).

If you had index groups A and B in index.ndx for ``g_sgangle``, you can use the
following command to compute the same angle with ``gmx gangle``::

  gmx gangle -n index.ndx -g1 vector/plane -group1 'group "A"' -g2 vector/plane -group2 'group "B"' -oav

You need to select either ``vector`` or ``plane`` for the ``-g1`` and ``-g2``
options depending on which one your index groups specify.

If you only had a single index group A in index.ndx and you used ``g_sgangle``
``-z`` or ``-one``, you can use::

  gmx gangle -n index.ndx -g1 vector/plane -group1 'group "A"' -g2 z/t0 -oav

For the distances, you can use :ref:`gmx distance` to compute one or more
distances as you want.  Both distances between centers of groups or individual
atoms are supported using the new selection syntax.

genbox
......

This tool has been split to :ref:`gmx solvate` and :ref:`gmx insert-molecules`.

tpbconv
.......

This tool has been renamed :ref:`gmx convert-tpr`.
