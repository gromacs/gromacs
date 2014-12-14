Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

energy split into four separate tools
"""""""""""""""""""""""""""""""""""""

The energy tool for analsys of energy files is split into four
separate
tools: one to print energies (called energy),
one to compute viscosity in a number of ways (called viscosity),
one to compute fluctuation properties such as compressibilty and
heat capacities (called fluctprops) and one to extract free energies
terms from energy files (called dhdl).
