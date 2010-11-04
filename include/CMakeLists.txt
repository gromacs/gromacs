# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${INCL_INSTALL_DIR}/gromacs
  COMPONENT development
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
