# Data: Nothing to build, just installation
install(DIRECTORY . DESTINATION ${DATA_INSTALL_DIR}
  COMPONENT data
  PATTERN "Makefile*" EXCLUDE
  PATTERN "CMake*" EXCLUDE
  PATTERN "cmake*" EXCLUDE
  PATTERN "template_doc.c" EXCLUDE
  PATTERN "Template.mak" EXCLUDE
  PATTERN "*~" EXCLUDE
)

install(FILES template/CMakeLists.txt.template
        DESTINATION ${DATA_INSTALL_DIR}
        RENAME template/CMakeLists.txt
        COMPONENT data)
