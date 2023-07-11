Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

mdrun now checks for excluded pairs beyond the cut-off with reaction-field and FEP
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With reaction-field electrostatics and free-energy calculations,
excluded atom pairs are not allowed to be beyond the Coulomb cut-off distance.
Now mdrun checks for this and throws an error when this occurs.

:issue:`4667`

enemat now prints correct headers when using ``-free`` or ``-eref`` options
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Fixed a long-standing bug when ``gmx enemat`` would output incorrect headers
to XVG.

:issue:`4812`
