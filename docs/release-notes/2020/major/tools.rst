Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

Fixed bug in gmx order -calcdist
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The reference position for the distance calculation was calculated
wrongly.

Improved grompp usability by rejecting more invalid .mdp lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Lines like

    ref-t 298
    = 0.1
    =

are now all rejected with a descriptive message, which will help
prevent some kinds of errors in constructing .mdp inputs. Note that an
.mdp parameter name with a missing value is still accepted, and leads
to the default behavior for that parameter.
