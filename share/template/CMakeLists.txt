install(FILES CMakeLists.txt.template
        DESTINATION ${DATA_INSTALL_DIR}/template
        RENAME CMakeLists.txt
        COMPONENT data)

install(FILES README template.c Makefile.pkg
        DESTINATION ${DATA_INSTALL_DIR}/template
        COMPONENT data)

install(FILES cmake/FindGROMACS.cmake
        DESTINATION ${DATA_INSTALL_DIR}/template/cmake
        COMPONENT data)

