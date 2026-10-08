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

PDB trajectory output now writes correct per-frame PBC type
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When writing PDB trajectories, the PBC type of the system was ignored
and instead guessed from the box. As a result, the ``CRYST1`` record
was written incorrectly for screw PBC (wrong space group and cell length).

Analysis tools now always process time values in double precision
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This avoids picking incorrect frames when times values are large.
