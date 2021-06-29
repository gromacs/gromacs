Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

grompp no longer modifies nstcomm
"""""""""""""""""""""""""""""""""

grompp will no longer set nstcomm, the interval for center of mass motion
removal, equal to nstcalcenergy when nstcomm < nstcalcenergy.
A note is still printed in that case.
