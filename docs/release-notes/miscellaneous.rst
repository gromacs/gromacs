Miscellaneous
^^^^^^^^^^^^^

grompp discourages use of constraints=all-bonds
===============================================
For force field parametrized with constraining only bonds involving
hydrogens, constraining all bonds should be avoided for correctness
and performance. grompp now issues a note when constraints=all-bonds
with Amber, CHARMM and OPLS-aa when time steps are smaller than 2.6 fs
and hydrogens are not replaced by virtual sites.
