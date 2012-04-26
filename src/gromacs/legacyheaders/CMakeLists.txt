# includes: Nothing to build, just configuration and installation
install(DIRECTORY . DESTINATION ${INCL_INSTALL_DIR}/gromacs/legacyheaders

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gmx_header_config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/gmx_header_config.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gmx_header_config.h
  DESTINATION ${INCL_INSTALL_DIR}/gromacs
  COMPONENT development
)

  COMPONENT development
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*.cmakein" EXCLUDE
)
