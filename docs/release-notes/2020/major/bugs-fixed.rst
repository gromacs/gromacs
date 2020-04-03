Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

gmx mdrun -append now requires that a checkpoint is found
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Previously ``gmx mdrun -append`` would start from the .tpr
configuration (and thus not append) when the checkpoint file was missing.
