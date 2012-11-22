add_executable(g_ana g_ana.cpp legacymodules.cpp)
target_link_libraries(g_ana libgromacs gmxana ${GMX_EXE_LINKER_FLAGS})
set_target_properties(g_ana PROPERTIES OUTPUT_NAME "g_ana${GMX_BINARY_SUFFIX}")
configure_file(CreateLinks.cmake.cmakein CreateLinks.cmake @ONLY)
set(CREATE_LINKS_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/CreateLinks.cmake)
# Based on documentation, CMake only supports symbolic links on Unix,
# although NTFS also has those (since Windows Vista; they require admin
# permissions to create).
if (UNIX)
    # TODO: Currently this rule is always run and recreates the symlinks;
    # it would be best to just adapt the rest of the build such that this is
    # never needed (for now, this makes the binary directory seem like it was
    # earlier).
    add_custom_target(links ALL
        ${CMAKE_COMMAND}
            -DBINARY_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            -DQUIETLY=TRUE
            -P ${CREATE_LINKS_SCRIPT}
        DEPENDS g_ana
        COMMENT "Creating symlinks for alternative binary names"
        VERBATIM)
endif ()
# TODO: Create man pages (currently not working since the old mechanism
# requires a target for each binary to which one can add a post-build role, and
# the symlinks do not have such targets).
# May be worth to wait until #685 and #969 are resolved; after these, it should
# be possible to handle everything uniformly instead of having several
# different CMake-side mechanisms.

install(TARGETS g_ana
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})
if (GMX_SYMLINK_OLD_BINARY_NAMES)
    install(SCRIPT ${CREATE_LINKS_SCRIPT} COMPONENT links)
endif ()
