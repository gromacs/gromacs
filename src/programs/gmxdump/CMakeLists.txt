file(GLOB GMXDUMP_SOURCES gmxdump.c ../main.cpp)

add_executable(gmxdump ${GMXDUMP_SOURCES})
gmx_add_man_page(gmxdump)
target_link_libraries(gmxdump ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
set_target_properties(gmxdump PROPERTIES OUTPUT_NAME "gmxdump${GMX_BINARY_SUFFIX}")
install(TARGETS gmxdump DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
