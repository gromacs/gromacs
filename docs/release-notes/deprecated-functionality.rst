.. _anticipated-changes:

Changes anticipated to |Gromacs| NEXT functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| NEXT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of virtual sites to replace aromatic rings in standard residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
These are thought to produce artefacts under some circumstances
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
Options such as ``-confout``, ``-resethway``, ``-resetstep`` are not
intended for use by regular mdrun users, so making them only available
with a dedicated tool is more clear. Also, this permits us to customize
defaults for e.g. writing files at the end of a simulation part in ways
that suit the respective mdrun and benchmark use cases, so ``-confout``
will no longer be required.

``gmx mdrun -nsteps``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The number of simulation steps described by the .tpr file can be
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
scheme. Deprecated in GROMACS 5.1
