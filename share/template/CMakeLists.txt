install(FILES CMakeLists.txt.template
        DESTINATION ${DATA_INSTALL_DIR}/template
        RENAME CMakeLists.txt
        COMPONENT development)

add_custom_command(OUTPUT gromacs
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/include gromacs)
# avoid removing of headers with make clean
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM TRUE)

add_executable(template template.c)
target_link_libraries(template gmx)
set_target_properties(template PROPERTIES INCLUDE_DIRECTORIES
  "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_BINARY_DIR}/src")
add_dependencies(template gromacs)

install(FILES README template.c Makefile.pkg
        DESTINATION ${DATA_INSTALL_DIR}/template
        COMPONENT development)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGROMACS.cmakein
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/FindGROMACS.cmake @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/FindGROMACS.cmake
        DESTINATION ${DATA_INSTALL_DIR}/template/cmake
        COMPONENT development)

