Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Dynamic pairlist generation for energy minimization
"""""""""""""""""""""""""""""""""""""""""""""""""""

With energy minimization, the pairlist, and domain decomposition when running
in parallel, is now performed when at least one atom has moved more than the
half the pairlist buffer size. The pairlist used to be constructed every step.
