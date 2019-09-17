.. _anticipated-changes:

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

Changes anticipated to |Gromacs| NEXT functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2020
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Support for 32bit architectures
"""""""""""""""""""""""""""""""
There are no current or planned large scale resources using 32bit architectures,
and we have no ability to properly test and evaluate them.


Functionality deprecated in |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of virtual sites to replace aromatic rings in standard residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
These are thought to produce artefacts under some circumstances
(unpublished results), were never well tested, are not widely used,
and we need to simplify pdb2gmx.

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
