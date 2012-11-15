add_executable(g_ana g_ana.cpp legacymodules.cpp)
target_link_libraries(g_ana libgromacs gmxana)
set_target_properties(g_ana PROPERTIES OUTPUT_NAME "g_ana${GMX_BINARY_SUFFIX}")
configure_file(CreateLinks.cmake.cmakein CreateLinks.cmake @ONLY)
set(CREATE_LINKS_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/CreateLinks.cmake)
if (UNIX)
    # TODO: Currently this rule is always run; it would be best to just adapt
    # the rest of the build such that this is never needed (for now, this makes
    # the binary directory seem like it was earlier).
    add_custom_target(links ALL
        ${CMAKE_COMMAND}
            -DBINARY_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DQUIETLY=TRUE
            -P ${CREATE_LINKS_SCRIPT}
        DEPENDS g_ana
        COMMENT "Creating symlinks for alternative binary names"
        VERBATIM)
    # TODO: Create man pages for legacy tools (if we want to have that working
    # before #969 is resolved, after which it should be possible to handle old
    # and new tools uniformly).
endif ()

install(TARGETS g_ana
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})
if (GMX_SYMLINK_OLD_BINARY_NAMES)
    install(SCRIPT ${CREATE_LINKS_SCRIPT} COMPONENT links)
endif ()
