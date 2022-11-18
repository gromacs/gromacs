Removed functionality
^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


The built-in viewer ``gmx view`` was removed
""""""""""""""""""""""""""""""""""""""""""""

There is little use and no tests of this functionality, and it was
deprecated in 2022.

:issue:`4296`

Some unmaintained utility scripts were removed
""""""""""""""""""""""""""""""""""""""""""""""

Several scripts in repository :file:`scripts/` directory were not installed with
the package, have not been maintained, and, as best we could tell, have not been
used in a long time.

:issue:`4639`
