file(GLOB GMXCHECK_SOURCES gmxcheck.c tpbcmp.c ../main.cpp)

add_executable(gmxcheck ${GMXCHECK_SOURCES})
gmx_add_man_page(gmxcheck)
target_link_libraries(gmxcheck ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
set_target_properties(gmxcheck PROPERTIES OUTPUT_NAME "gmxcheck${GMX_BINARY_SUFFIX}")
install(TARGETS gmxcheck DESTINATION ${BIN_INSTALL_DIR} COMPONENT runtime)
