# Data: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${DATA_INSTALL_DIR}
  PATTERN "Makefile*" EXCLUDE
  PATTERN ".svn" EXCLUDE
  PATTERN "CVS" EXCLUDE
  PATTERN ".cvsignore" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
