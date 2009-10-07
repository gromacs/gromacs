# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION include/gromacs/thread_mpi/atomic
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
