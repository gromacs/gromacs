include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

file(GLOB PDB2GMX_SOURCES *.c)

add_executable(pdb2gmx ${PDB2GMX_SOURCES})
gmx_add_man_page(pdb2gmx)
target_link_libraries(pdb2gmx ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
set_target_properties(pdb2gmx PROPERTIES OUTPUT_NAME "pdb2gmx${GMX_BINARY_SUFFIX}")
install(TARGETS pdb2gmx DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
