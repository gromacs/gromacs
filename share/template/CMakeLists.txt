install(FILES CMakeLists.txt.template
        DESTINATION ${DATA_INSTALL_DIR}/template
        RENAME CMakeLists.txt
        COMPONENT development)

install(FILES README template.c Makefile.pkg
        DESTINATION ${DATA_INSTALL_DIR}/template
        COMPONENT development)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGROMACS.cmakein
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/FindGROMACS.cmake @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/FindGROMACS.cmake
        DESTINATION ${DATA_INSTALL_DIR}/template/cmake
        COMPONENT development)

