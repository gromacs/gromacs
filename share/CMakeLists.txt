# Data: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${DATA_INSTALL_DIR}
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
