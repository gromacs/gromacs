Removed functionality
^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

gmx anadock
"""""""""""
The gmx anadock tool was removed since it does not belong in gromacs
(it analyzes AutoDock outputs).

gmx dyndom
""""""""""
The gmx dyndom tool was removed since it does not belong in gromacs
(it analyzes DynDom outputs).

gmx morph
"""""""""
The gmx morph tool was removed since it yields non-physical structures
that can easily be done by a script.
