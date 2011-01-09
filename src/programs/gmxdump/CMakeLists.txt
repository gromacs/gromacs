include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

file(GLOB GMXDUMP_SOURCES gmxdump.c)

add_executable(gmxdump ${GMXDUMP_SOURCES})
target_link_libraries(gmxdump ${GMX_EXTRA_LIBRARIES} libgromacs)
set_target_properties(gmxdump PROPERTIES OUTPUT_NAME "gmxdump${GMX_BINARY_SUFFIX}")
install(TARGETS gmxdump DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
