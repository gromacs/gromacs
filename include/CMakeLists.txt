# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION include/gromacs
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
