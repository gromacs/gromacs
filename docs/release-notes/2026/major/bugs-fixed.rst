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
