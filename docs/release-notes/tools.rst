Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added option -water tips3p to pdb2gmx.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2272`

Removed incorrect comment for CHARMM tips3p
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Removed CHARMM tips3p performance warning in :ref:`gmx pdb2gmx` input file,
since the performance loss is negligible with the cutoff-scheme=Verlet.

Fixed :ref:`gmx check` for tprs with different numbers of atoms
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2279`

Split off the NMR related analyses from :ref:`gmx energy`.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A new tool :ref:`gmx nmr` is created by straight copying code from
:ref:`gmx energy` to a new tool. The reason is to reduce complexity.

A few cleanups are introduced to pass the valgrind memory
test.

Added references the :ref:`gmx nmr` in the manual.

Avoided :ref:`gmx grompp` charge warning from merely rounding error
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Even though the :ref:`gmx grompp` total charge check uses double for summation,
there are rounding errors for each charge when charges are stored
in single precision. Now the charge check rounds the net charge of
molecules to integer when the difference is less than the maximum
possible sum of charge rounding errors.

Fixes :issue:`2192`

Made duplicate atoms in bondeds an error in :ref:`gmx grompp`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Having duplicate atom indices in bonded interactions used to be only
a warning. But since in nearly all cases this will lead to issues,
this is now a error, except for angle restraints where it can be
useful so there it is now a note.

:issue:`2141`

Clarified :ref:`gmx editconf` help text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
It is possible that users can confuse -c with -center so this
patch makes it clear that -center doesn't do anything unless the
user really wants to shift the center of the system away from the
middle of the box.

Fixes :issue:`2171`

Decreased memory usage in :ref:`gmx traj` and :ref:`gmx trjconv`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Made :ref:`gmx grompp` -r obligatory with position restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With position restraints it would often occur that users accidentally
used equilibrated coordinates instead of the original coordinates for
position restraint coordinates due to :ref:`gmx grompp` -r defaulting
to -c. Now -r always need to be supplied with position restraints,
but using the same file name as with -c will reproduce the old
behavior.

Added selection-enabled :ref:`gmx traj`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
For now, this tool only plots coordinates, velocities, and forces for
selections, so it should provide a full replacement for -ox, -ov, -of,
-com, and -mol from :ref:`gmx traj`.
