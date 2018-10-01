Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Split off the NMR related analyses from :ref:`gmx energy`.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A new tool :ref:`gmx nmr` is created by straight copying code from
:ref:`gmx energy` to a new tool. The reason is to reduce complexity.

A few cleanups are introduced to pass the valgrind memory test.

Added references the :ref:`gmx nmr` in the manual.

Added selection-enabled :ref:`gmx trajectory`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
For now, this tool only plots coordinates, velocities, and forces for
selections, so it should provide a full replacement for -ox, -ov, -of,
-com, and -mol from :ref:`gmx traj`.

Decreased memory usage in :ref:`gmx traj` and :ref:`gmx trjconv`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Made TNG writing work with multiple identical steps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Introduce a wrapper structure around TNG so we detect and correct for
cases when writing multiple frames with the same step, or non-zero
initial steps to TNG files.  This will avoid frames overwriting each
other, and make sure the time per frame is correct.

:issue:`2189`

Improved frame time/step handling in :ref:`gmx trjconv`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Store the exact step in PDB/GRO file headers, and be more careful
about not claiming to have time or step information when it was not
available.  This change will avoid some of the problems described in
:issue:`2189`, but it does not yet properly fix the issue in the TNG
library.

:issue:`2189`

Fixed :ref:`gmx trjconv` to always dump at correct time
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Set frame timestep before starting the loop by reading first two
frames and rewinding, and make sure we always write something to the
dump output based on best-guess (if there is at least one input frame
present).

:issue:`1832`

Clarified :ref:`gmx editconf` help text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
It is possible that users can confuse -c with -center so this
patch makes it clear that -center doesn't do anything unless the
user really wants to shift the center of the system away from the
middle of the box.

Fixes :issue:`2171`

Added option -water tips3p to pdb2gmx.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2272`

Removed incorrect comment for CHARMM tips3p
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Removed CHARMM tips3p performance warning in :ref:`gmx pdb2gmx` input file,
since the performance loss is negligible with the
:mdp-value:`cutoff-scheme=Verlet`.

Avoided :ref:`gmx grompp` charge warning from merely rounding error
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Even though the :ref:`gmx grompp` total charge check uses double for summation,
there are rounding errors for each charge when charges are stored
in single precision. Now the charge check rounds the net charge of
molecules to integer when the difference is less than the maximum
possible sum of charge rounding errors.

Fixes :issue:`2192`

Improved pdb2gmx for nonstandard residue types
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If explicit non-blank chain identifiers are set, it will now be a hard
error if the residue types in each chain do not match. For blank chain
ID we still need to allow detection of non-chain parts, but this case
too now provides more explicit output information.

:issue:`2370`

Allowed empty lines in hdb files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Skip lines that consist only of whitespace. Not a universal solution
for fixing hdb files, but better than the user getting very strange
error messages that don't say anything about whitespace.

:issue:`2028`

Changed to no longer require matching names between rtp and tdb files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This was only documented in the source. It's a remnant from the days
when all force fields were in the same directory, and no longer
necessary. With this change we will properly match all termini to all
amino acids.

:issue:`2026`
:issue:`2027`

Made duplicate atoms in bondeds an error in :ref:`gmx grompp`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Having duplicate atom indices in bonded interactions used to be only
a warning. But since in nearly all cases this will lead to issues,
this is now a error, except for angle restraints where it can be
useful so there it is now a note.

:issue:`2141`

Made :ref:`gmx grompp` -r obligatory with position restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With position restraints it would often occur that users accidentally
used equilibrated coordinates instead of the original coordinates for
position restraint coordinates due to -r defaulting
to -c. Now -r always need to be supplied with position restraints,
but using the same file name as with -c will reproduce the old
behavior.

Fixed :ref:`gmx msd` when using COM removal and molecules
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Changed order of code to actually assign correct coordinates before
copying the data, and modified data structure size when using COM
removal and individual molecules.

:issue:`2043`

Fixed index error in :ref:`gmx chi`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
An error in the index construction could lead to segfaults. However,
the actual indices were correct, so it should not have produced any
incorrect results.

:issue:`1814`

Fixed :ref:`gmx grompp` complexity for large exclusion orders
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
To avoid exploding computational complexity for highly connected
molecules with large values for excluded neighbors, avoid adding a
neighbor to the temporary nnb structure if it is already present as a
lower-order neighbor.

:issue:`2260`

Fixed :ref:`gmx density` for non-mass calculations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now always use mass and never charge/electron density to center
systems.

:issue:`2230`

Fixed :ref:`gmx check` for tprs with different numbers of atoms
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2279`
