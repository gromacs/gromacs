Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Renamed all NBNXN environment variables to NBNXM
""""""""""""""""""""""""""""""""""""""""""""""""

All environment variables starting with ``GMX_NBNXN`` have been renamed to start with ``GMX_NBNXM`` for consistency with the algorithm's name.

Replaced usage of custom Bohr radius value in gmx spatial with the common value from units
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :ref:`gmx spatial` command used to have its own definition of Bohr radius. For consistency
with other parts of |Gromacs|, it now uses the definition of Bohr radius from the same source as
the rest of the code. Notably, the value of the constant changed from ``0.529177249`` (IUPAC 1999)
to ``0.529177210903`` (NIST 2018).

AMBER19SB and AMBER14SB force fields now use IUPAC standard hydrogen names
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:ref:`gmx pdb2gmx` now recognizes the hydrogen names from the force field rather than always
starting from 1. This means that the AMBER19SB and AMBER14SB force fields now use the original
hydrogen names, (i.e. HB2 and HB3 for methylene hydrogens instead of HB1 and HB2). This is more
consistent with the naming in the original Amber force field files and with the IUPAC standard
for hydrogen names.