add_subdirectory(template)

set(HTML_DIR ${CMAKE_CURRENT_BINARY_DIR}/html/online)
file(MAKE_DIRECTORY ${HTML_DIR})
add_custom_target(html
    gmx -quiet help -export html
    WORKING_DIRECTORY ${HTML_DIR}
    COMMENT "Generation HTML help for programs"
    VERBATIM)

install(DIRECTORY . DESTINATION ${DATA_INSTALL_DIR}
  COMPONENT data
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "*~" EXCLUDE
  PATTERN "template" EXCLUDE
)
