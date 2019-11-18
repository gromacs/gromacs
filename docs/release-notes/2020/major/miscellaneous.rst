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

Provide checksum to validate release tarballs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Released versions of |Gromacs| will now provide a checksum calculated
from the files participating in building the binaries. When building
|Gromacs| from the tarball, the files will be checksummed again and
compared against the checksum generated during the release build. If the
checksums don't match, the version string is modified to indicate that
the source tree has been modified, and the information is printed in the
log files for the users. If checksumming has not been possible (either due
to missing Python during installation, or because the original checksum file
is missing), this is indicated through a different version string.

:issue:`2128`

