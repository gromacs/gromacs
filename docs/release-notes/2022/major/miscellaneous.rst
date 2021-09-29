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

Bonded atom types names can now start with a digit
""""""""""""""""""""""""""""""""""""""""""""""""""

Bonded atom types names in topologies were not allowed to start with a number.
Now all names are supported that contain at least one non-digit character.

:issue:`4120`
