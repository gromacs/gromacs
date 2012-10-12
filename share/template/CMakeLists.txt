add_executable(template template.cpp)
target_link_libraries(template libgromacs)
install(FILES CMakeLists.txt.template
        DESTINATION ${DATA_INSTALL_DIR}/template
        RENAME CMakeLists.txt
        COMPONENT development)

install(FILES README template.cpp Makefile.pkg
        DESTINATION ${DATA_INSTALL_DIR}/template
        COMPONENT development)

install(FILES cmake/FindGROMACS.cmake
        DESTINATION ${DATA_INSTALL_DIR}/template/cmake
        COMPONENT development)

