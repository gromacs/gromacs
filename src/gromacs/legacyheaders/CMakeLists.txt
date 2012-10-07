# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${INCL_INSTALL_DIR}/gromacs/legacyheaders
  COMPONENT development
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*.cmakein" EXCLUDE
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/version.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h
  DESTINATION ${INCL_INSTALL_DIR}/gromacs
  COMPONENT development
)
