# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${INCL_INSTALL_DIR}/gromacs
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
