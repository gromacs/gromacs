add_executable(g_options options.c)
target_link_libraries(g_options ${GMX_EXTRA_LIBRARIES} libgromacs)
set_target_properties(g_options PROPERTIES OUTPUT_NAME "g_options${GMX_BINARY_SUFFIX}")
