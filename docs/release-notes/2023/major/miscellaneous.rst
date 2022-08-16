Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Fix documentation issues for restricted bending potential 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The actual implementation in the code was correct, but the
manual section for the restricted bending potential had an
extra factor 2 for the force inherited from the 2013
Bulacu JCTC paper, and the journals for two of the references
had been swapped. No changes to any simulation results.

:issue:`4568`
