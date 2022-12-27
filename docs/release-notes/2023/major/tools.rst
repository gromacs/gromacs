Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


``gmx do_dssp`` replaced by native implementation of DSSP algorithm
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``gmx do_dssp`` replaced by native implementation of DSSP algorithm,
version 4. Results for version 4 was formerly available with ``gmx
do_dssp -ver 4``. The new tool is called ``gmx dssp``.

