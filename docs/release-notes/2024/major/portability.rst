Portability
^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Always use the Boost version bundled with |Gromacs|
"""""""""""""""""""""""""""""""""""""""""""""""""""

Boost 1.83 is known to have compatibility issues when using Clang
compiler on FreeBSD and Linux. This changes to always use the bundled
Boost version, even if another version is present on the system.

:issue:`4893`
