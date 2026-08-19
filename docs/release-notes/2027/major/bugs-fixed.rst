Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

PDB with a box of size 1 Angstrom are interpreted as having no PBC
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

PDB files with structures not determined from X-ray should have a unit-cell with P=1
and dimensions of 1 Angstrom to indicate no periodic boundary conditions.
|Gromacs| now interprets such structures as not having PBC.

:issue:`4645`, :issue:`5679`
