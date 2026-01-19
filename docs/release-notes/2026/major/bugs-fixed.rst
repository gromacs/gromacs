Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


Collective variables module (Colvars) update
""""""""""""""""""""""""""""""""""""""""""""

The (`Colvars <https://colvars.github.io>`_) library for enhanced sampling simulations included
in |Gromacs| has been updated to version 2025-10-13.

This update brings many fixes, including some errors in gradient calculations, mostly affecting
rarely used cases (see `this page for details <https://github.com/Colvars/colvars/wiki/Gradient-evaluation-bugs-fixed-in-September-2025>`_).
A complete list of changes can be found `here <https://gitlab.com/gromacs/gromacs/-/merge_requests/5397>`_.

Tpr-file versioning fixed
"""""""""""""""""""""""""

Older installed |Gromacs| versions are intended to be able to do at least basic reading of
.tpr files written by newer versions of |Gromacs|, unless made impossible by changes to the
.tpr format. This happened in |Gromacs| 2025, but was not accompanied an appropriate increment
of the version numbers, resulting in incorrect and misleading error messages.
This cannot be fixed retroactively for installed |Gromacs| versions and extant .tpr files,
but at least the problem conditions will no longer continue to arise for newly created
.tpr files.

:issue:`5334`

Fixed missing non-bonded interactions with direct halo communication
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When using the experimental direct halo communication feature combined with 8-wide SIMD
non-bonded kernels and OpenMP threading, non-bonded interactions could be missing.

:issue:`5509`

Allow atoms involved intermolecular-exclusion to be perturbed
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:issue:`5527`

Added check for constructing atoms of virtual sites
"""""""""""""""""""""""""""""""""""""""""""""""""""

Constructing atoms for virtual sites can themselves be virtual sites, but only when
those constructing atoms are virtual sites are higher up in the function type list
(i.e. simpler constructions). This was documented in the manual. Now``grompp``
and ``mdrun`` will throw an error when these restrictions are violated.

:issue:`5535`
