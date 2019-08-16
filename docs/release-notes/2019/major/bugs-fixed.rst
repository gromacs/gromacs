Bugs fixed
^^^^^^^^^^

Fix type bug in trilinic DD code
""""""""""""""""""""""""""""""""""""""""""""""""""

Fix bug with unusual off-diagonal elements communicating too few atoms.

Ensure domains are large enough for atom motion
""""""""""""""""""""""""""""""""""""""""""""""""""

Domain decomposition now makes sure that domains will always be large
enough so that atoms will not move across additional domains.

:issue:`2614`

Velocity Verlet integrators output energy averages from correct steps
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Velocity Verlet integrators would accumulate energies for writing
averages to the energy file when step-1 was a multiple of nstcalcenergy.
This has now been corrected to step being a multiple of nstcalcenergy.
Note that although this (slightly) changes the reported averages,
the averages were not incorrect.

:issue:`2718`

Fix chainsep behaviour of pdb2gmx
""""""""""""""""""""""""""""""""""""""""""""""""""

:issue:`2577`

grompp correctly checks nstexpanded against nstcalcenergy
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With expanded ensemble, but without free-energy perturbation, grompp
would not check if nstexpanded was a multiple of nstcalcenergy.
If the latter was not the case, results might have been incorrect.

:issue:`2714`

Issue with do_dssp and unknown residues
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :ref:`do_dssp <gmx do_dssp>` tool would fail with unknown residues,
as well as have issues on Windows.

:issue:`2599`

