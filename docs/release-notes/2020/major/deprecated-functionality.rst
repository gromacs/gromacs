.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

The core |Gromacs| team wants to let users and downstream developers
know about impending changes so that disruption is minimized. Do get
in touch if you feel something inappropriate is planned!

Deprecated functionality often remains in |Gromacs| for a year or
more, but this should not be relied upon.

Changes anticipated to |Gromacs| 2020 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``gmx mdrun -membed``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The feature for embedding a protein in a membrane will be retained,
but probably in a different form, such as ``gmx membed``.

``gmx mdrun -rerun``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The feature for computing potential energy quantities from a
trajectory will be retained, but probably in a different form, such as
``gmx rerun`` and ``gmx test-particle-insertion``.

Integrator .mdp options will only contain dynamical integrators
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Energy minimization will be accessed in a differt form, perhaps with
``gmx minimize`` and interpret an .mdp field for which minimizer to
use. Normal-mode analysis may be accessed with e.g. ``gmx
normal-modes``. The command-line help for these tools will then
be better able to document which functionality is supported when.

Much functionality in ``trjconv``, ``editconf``, ``eneconv`` and ``trjcat``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The functionality in such tools is being separated to make it
available in composable modules, that we plan to make available as
simpler tools, and eventually via the |Gromacs| API that is under
development.

``gmx do_dssp`` to be replaced
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This tool is deprecated, because it is problematic for some users to
obtain and install a separate DSSP binary, so we plan to replace the
implementation at some point with a native implementation, likely
based upon xssp, and make it available under a new gmx tool name.

Functionality deprecated in |Gromacs| 2020
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Support for 32bit architectures
"""""""""""""""""""""""""""""""
:issue:`3252` There are no current or planned large scale resources using 32bit architectures,
and we have no ability to properly test and evaluate them.

Free-energy soft-core power 48
""""""""""""""""""""""""""""""
:issue:`3253` Free-energy soft-core power 48 is almost never used and is therefore deprecated.

Support for Armv7
"""""""""""""""""
:issue:`2990` There are several issues with current code for the architecture, and we don't
have the resources for support and fix issues related to it. As the architecture has no
large HPC impact it is thus deprecated.

Functionality deprecated in |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of virtual sites to replace aromatic rings in standard residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3254` These are thought to produce artefacts under some circumstances
(unpublished results), were never well tested, are not widely used,
and we need to simplify pdb2gmx.

Benchmarking options only available with ``gmx benchmark``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3255` Options such as ``-confout``, ``-resethway``, ``-resetstep`` are not
intended for use by regular mdrun users, so making them only available
with a dedicated tool is more clear. Also, this permits us to customize
defaults for e.g. writing files at the end of a simulation part in ways
that suit the respective mdrun and benchmark use cases, so ``-confout``
will no longer be required.

``gmx mdrun -nsteps``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3256` The number of simulation steps described by the .tpr file can be
changed with ``gmx convert-tpr``, or altered in .mdp file before the
call to ``gmx grompp``. The convenience of this mdrun option was
outweighted by the doubtful quality of its implementation, no clear
record in the log file, and lack of maintenance.
