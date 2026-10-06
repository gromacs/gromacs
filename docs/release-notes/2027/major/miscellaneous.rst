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

Velocities of special or frozen particles handled differently
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Formerly :ref:`gmx grompp` zeroed velocities of shells and virtual
sites when generating velocities, and then later :ref:`gmx mdrun`
zeroed them again regardless of their origin, and also zeroed velocity
components of particles with frozen dimensions. Now :ref:`gmx mdrun`
does nothing, and :ref:`gmx grompp` zeroes frozen velocity
components. Users will now see that velocities of special particles
passed to :ref:`gmx grompp` are retained by :ref:`gmx mdrun`.

:issue:`5714`

AMBER19SB and AMBER14SB force fields now use IUPAC standard hydrogen names
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:ref:`gmx pdb2gmx` now recognizes the hydrogen names from the force field rather than always
starting from 1. This means that the AMBER19SB and AMBER14SB force fields now use the original
hydrogen names, (i.e. HB2 and HB3 for methylene hydrogens instead of HB1 and HB2). This is more
consistent with the naming in the original Amber force field files and with the IUPAC standard
for hydrogen names.

Warnings previously changing user inputs converted to fatal errors
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Various warnings in analysis tools previously would warn the user that input combinations were
incompatible and the desired output would not be produced. To prevent confusion or the need to
parse lengthy log files in order to verify the outputs of the tools, these warnings were converted
to fatal errors so that only valid combinations of inputs will result in successful execution.

:issue:`5626`
