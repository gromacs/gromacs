# includes: Nothing to build, just configuration and installation
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gmxconfig.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/gmxconfig.h)

install(DIRECTORY . DESTINATION ${INCL_INSTALL_DIR}/gromacs
  COMPONENT development
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
)
