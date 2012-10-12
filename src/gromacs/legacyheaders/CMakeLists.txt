# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${INCL_INSTALL_DIR}/gromacs/legacyheaders
  COMPONENT development
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*.cmakein" EXCLUDE
)
