Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Fixed exported libgromacs CMake target
""""""""""""""""""""""""""""""""""""""

Update the exported libgromacs CMake target to not depend on non-
existing include paths and add GMX_DOUBLE define to interface
definitions. The target now gets exported into the Gromacs namespace.

:issue:`3468`
