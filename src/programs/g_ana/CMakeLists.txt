add_executable(g_ana g_ana.cpp)
target_link_libraries(g_ana libgromacs)
set_target_properties(g_ana PROPERTIES OUTPUT_NAME "g_ana${GMX_BINARY_SUFFIX}")

install(TARGETS g_ana
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})
if (GMX_SYMLINK_OLD_BINARY_NAMES)
    configure_file(CreateLinks.cmake.cmakein CreateLinks.cmake @ONLY)
    install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/CreateLinks.cmake COMPONENT links)
endif ()
