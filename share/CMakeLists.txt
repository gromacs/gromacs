# Data: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${DATA_INSTALL_DIR}
  COMPONENT data
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
  PATTERN "template" EXCLUDE
)

add_subdirectory(template)
