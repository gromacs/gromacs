Deprecated functionality
------------------------

Changes anticipated to |Gromacs| 2019 functionality
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
normal-modes``. The command-line help for these tools will thenx
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

Functionality deprecated in |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of virtual sites to replace aromatic rings in standard residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3254` These are thought to produce artefacts under some circumstances
(unpublished results), were never well tested, are not widely used,
and we need to simplify pdb2gmx.

``gmx mdrun -gcom``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This feature sometimes overrides the effects of various .mdp settings
in a way that is difficult to understand and report. A user who wants
to do communication between PP ranks less often should choose their
``nst*`` mdp options accordingly.

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

Functionality deprecated before |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This functionality has been declared as deprecated in previous versions
of |Gromacs|, but has not yet been removed.

The group cutoff scheme
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
All remaining aspects of the group cutoff scheme will be removed, once
a few remaining features (e.g. tabulated interactions, energy-group
exclusions, and vacuum simulations) are available with the Verlet
scheme. Deprecated in |Gromacs| 5.1

QM/MM support for ORCA, GAMESS and MOPAC
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
These interfaces are untested, and no maintainer has been found for them.
Deprecated in |Gromacs| 2018.
