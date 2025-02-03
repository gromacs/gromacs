Removed functionality
^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Removed support for planetary simulations
"""""""""""""""""""""""""""""""""""""""""

``GMX_DO_GALACTIC_DYNAMICS`` setting is removed.
Setting a negative :mdp:`epsilon-r` value now always raises an error.

:issue:`5223`

