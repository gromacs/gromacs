# Man pages: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${MAN_INSTALL_DIR}
  PATTERN "Makefile*" EXCLUDE
  PATTERN ".svn" EXCLUDE
  PATTERN "CVS" EXCLUDE
  PATTERN ".cvsignore" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
