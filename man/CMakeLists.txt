# Man pages: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${MAN_INSTALL_DIR}
  COMPONENT data
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
