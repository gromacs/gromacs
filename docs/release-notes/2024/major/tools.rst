Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

XVG output from ``gmx rdf`` now uses 6 decimal places
"""""""""""""""""""""""""""""""""""""""""""""""""""""

The output from ``gmx rdf`` now uses more decimal places in order to
avoid rounding issues. These issues led to perceived erroneous shifts in
the results.

:issue:`4647`
