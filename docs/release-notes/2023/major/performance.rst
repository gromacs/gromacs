Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Increased default T- and P-coupling intervals
"""""""""""""""""""""""""""""""""""""""""""""

The default maximum values temperature and pressure coupling intervals
have been increased from 10 to 100 steps. These values are used when
the default value of -1 is specified in the mdp file and a lower value
is used when required for accurate integration. The improves the performance
of both GPU runs and parallel runs.
