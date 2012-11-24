include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

file(GLOB X2TOP_SOURCES g_x2top.c nm2type.c)

add_executable(g_x2top ${X2TOP_SOURCES})
gmx_add_man_page(g_x2top)
target_link_libraries(g_x2top ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
set_target_properties(g_x2top PROPERTIES OUTPUT_NAME "g_x2top${GMX_BINARY_SUFFIX}")
install(TARGETS g_x2top DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
