include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

set(exename g_x2top)
configure_file(../main.cpp main.cpp)
file(GLOB X2TOP_SOURCES g_x2top.c nm2type.c ${CMAKE_CURRENT_BINARY_DIR}/main.cpp)

add_executable(g_x2top ${X2TOP_SOURCES})
gmx_add_man_page(g_x2top)
target_link_libraries(g_x2top ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
set_target_properties(g_x2top PROPERTIES OUTPUT_NAME "g_x2top${GMX_BINARY_SUFFIX}")
install(TARGETS g_x2top DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
