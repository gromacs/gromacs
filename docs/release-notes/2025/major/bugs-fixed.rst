Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Fix double counting of dispersion correction with test particle insertion
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The dispersion correction of the test particle with the system was double
counted.

:issue:`5231`

``gmx grompp`` now allows test particle insertion of an identical molecule
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``gmx grompp`` merges molecule blocks in the topology with the same molecule
type. With TPI this prevented insertion of a molecule into a liquid of the
same molecules.

