# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION include/gromacs
  PATTERN "Makefile*" EXCLUDE
  PATTERN ".svn" EXCLUDE
  PATTERN "CVS" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN ".cvsignore" EXCLUDE
  PATTERN "*~" EXCLUDE
)
