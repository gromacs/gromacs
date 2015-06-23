Command-line reference
======================

.. toctree::
   :hidden:
   :glob:

   /onlinehelp/gmx
   /onlinehelp/gmx-*

|Gromacs| includes many tools for preparing, running and analysing
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

.. toctree::

   /onlinehelp/selections

.. _command-changes:

Command changes between versions
--------------------------------

Starting from |Gromacs| 5.0, some of the analysis commands (and a few other
commands as well) have changed significantly.
In the process, some old analysis tools have been removed in favor of more
powerful functionality that is available through an alternative tool.
This page documents how to perform different tasks that were possible with the
old tools with the new set of tools.

Many of the new tools mentioned below now accept selections through one or more
command-line options instead of prompting for an index group.  Please see
:doc:`/onlinehelp/selections` additional information on how to use the
selections.

5.0
^^^

General
.......

Version 5.0 introduced the :command:`gmx` wrapper binary.
For backwards compatibility, this version still creates symbolic links by default for
old tools: e.g., ``g_order <options>`` is equivalent to ``gmx order <options>``, and
``g_order`` is simply a symbolic link on the file system.

g_bond
......

This tool has been removed in 5.0. A replacement is :ref:`gmx distance`.

You can provide your existing index file to :ref:`gmx distance`, and it will
calculate the same distances.  The only difference is that ``-blen`` and
``-tol`` options have different default values, and that you can control the
output histogram with ``-binw``.  ``-aver`` and ``-averdist`` options are not
present.  Instead, you can choose between the different things to calculate
using ``-oav`` (corresponds to ``-d`` with ``-averdist``), ``-oall``
(corresponds to ``-d`` without ``-averdist``), ``-oh`` (corresponds to ``-o``
with ``-aver``), and ``-oallstat`` (corresponds to ``-l`` without ``-aver``).
You can produce any combination of output files.  Compared to ``g_bond``,
``gmx distance -oall`` is currently missing labels for the output columns.

g_dist
......

This tool has been removed in 5.0.  A replacement is :ref:`gmx distance` (for
most options) or :ref:`gmx select` (for ``-dist`` or ``-lt``).

If you had index groups A and B in :file:`index.ndx` for ``g_dist``, you can use the
following command to compute the same distance with gmx distance::

  gmx distance -n index.ndx -select 'com of group "A" plus com of group "B"' -oxyz -oall

The ``-intra`` switch is replaced with ``-nopbc``.

If you used ``-dist D``, you can do the same calculation with ``gmx select``::

  gmx select -n index.ndx -select 'group "B" and within D of com of group "A"' -on/-oi/-os/-olt

You can select the output option that best suits your post-processing needs
(``-olt`` is a replacement for ``g_dist -dist -lt``)

g_sas
.....

This tool has been rewritten in 5.0, and renamed to :ref:`gmx sasa` (the
underlying surface area calculation algorithm is still the same).

The main difference in the new tool is support for selections.  Instead of
prompting for an index group, a (potentially dynamic) selection for the
calculation can be given with ``-surface``.  Any number of output groups can be
given with ``-output``. The total area of the ``-surface`` group is now always
calculated. Please see ``gmx sasa -h``.

The tool no longer automatically divides the surface into hydrophobic and
hydrophilic areas, and there is no ``-f_index`` option.  The same effects can
be obtained by defining suitable selections for ``-output``.  If you want
output that contains the same numbers as with the old tool for a calculation
group A and output group B, you can use ::

  gmx sasa -surface 'group "A"' -output '"Hydrophobic" group "A" and charge {-0.2 to 0.2}; "Hydrophilic" group "B" and not charge {-0.2 to 0.2}; "Total" group "B"'

Solvation free energy estimates are now calculated only if separately requested
with ``-odg``, and are written into a separate file.

Output option ``-i`` for a position restraint file is not currently implemented
in the new tool, but would not be very difficult to add if requested.

g_sgangle
.........

This tool has been removed in 5.0. A replacement is :ref:`gmx gangle` (for angle
calculation) and :ref:`gmx distance` (for ``-od``, ``-od1``, ``-od2``).

If you had index groups A and B in index.ndx for ``g_sgangle``, you can use the
following command to compute the same angle with ``gmx gangle``::

  gmx gangle -n index.ndx -g1 vector/plane -group1 'group "A"' -g2 vector/plane -group2 'group "B"' -oav

You need to select either ``vector`` or ``plane`` for the ``-g1`` and ``-g2``
options depending on which one your index groups specify.

If you only had a single index group A in index.ndx and you used ``g_sgangle``
``-z`` or ``-one``, you can use::

  gmx gangle -n index.ndx -g1 vector/plane -group1 'group "A"' -g2 z/t0 -oav

For the distances, you can use :ref:`gmx distance` to compute one or more
distances as you want. Both distances between centers of groups or individual
atoms are supported using the new selection syntax.

genbox
......

This tool has been split to :ref:`gmx solvate` and :ref:`gmx insert-molecules`.

tpbconv
.......

This tool has been renamed :ref:`gmx convert-tpr`.
