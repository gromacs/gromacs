# Man pages: Nothing to build, just installation
# Nothing to do for sources obtained from Git (man pages are auto-build (see 
# cmake/BuildManPages.cmake). But source packages contain man pages in this
# folder and those need to be installed here.
install(DIRECTORY . DESTINATION ${MAN_INSTALL_DIR}
  COMPONENT data
  FILES_MATCHING PATTERN "*.[1-9]"
  PATTERN "CMakeFiles" EXCLUDE
)
