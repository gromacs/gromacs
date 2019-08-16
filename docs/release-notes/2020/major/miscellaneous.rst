Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

grompp now warns if macros in mdp "define" field are unused in topology
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Macros defined in the mdp (with e.g. -DPOSRES) now cause a warning
in grompp if they are not encountered while parsing the topology file

:issue:`1975`

Introduced CMake variable GMX_VERSION_STRING_OF_FORK
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
To help users and developers understand which version of |Gromacs| is
being used, anybody providing a forked version of |Gromacs| shall set 
GMX_VERSION_STRING_OF_FORK in the source code (or if necessary when 
running CMake). It will then appear in the log file and users will know
which version and fork of the code produced the result.
