Removed functionality
^^^^^^^^^^^^^^^^^^^^^

Removed constant-acceleration groups support
""""""""""""""""""""""""""""""""""""""""""""
This code has been broken since before GROMACS 4.6, so it has been
removed.

:issue:`1354`

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Removed mdrun-only build configuration
""""""""""""""""""""""""""""""""""""""

The need for the mdrun-only build of |Gromacs| has expired, as it has
the same set of dependencies as regular |Gromacs|. It was deprecated
in GROMACS 2021. Removing it will simplify maintenance, testing,
documentation, installation, and teaching new users.

:issue:`3808`

Removed support for x86 MIC platform
""""""""""""""""""""""""""""""""""""

This platform is nearly dead and is no longer supported. The KNL
platform is unaffected by this change.

:issue:`3891`
