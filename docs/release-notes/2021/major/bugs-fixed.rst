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

Fixed unsolicited changing of atom names in pdb file
""""""""""""""""""""""""""""""""""""""""""""""""""""

Remove functions to change atoms names when reading 
and writing pdb files. This affected naming of
H atoms in particular.

:issue:`3469`

Corrected AWH initial histogram size
""""""""""""""""""""""""""""""""""""

The initial histogram size for AWH biases depended (weakly) on the force
constant. This dependence has been removed, which increases the histogram
size by a about a factor of 3. In practice this has only a minor effect
on the time to solution. For multiple dimensions, the histogram size was
underestimated, in particular with a combination of slower and faster
dimensions. The, now simplified, formula for the initial histogram size is
given in the reference manual.

:issue:`3751`
