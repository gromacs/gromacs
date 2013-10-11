configure_file(config.h.cmakein config.h)
configure_file(buildinfo.h.cmakein buildinfo.h)

if (BUILD_TESTING)
    add_custom_target(tests)
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(external/gmock-1.6.0)
    endif (GMX_BUILD_UNITTESTS)
    include(testutils/TestMacros.cmake)
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(testutils)
    endif (GMX_BUILD_UNITTESTS)
endif (BUILD_TESTING)

add_subdirectory(gromacs)
add_subdirectory(programs)

if(NOT GMX_FAHCORE)
  add_subdirectory(contrib)
endif()

# Set global variable from the local variable's value
set(TARGETS_NEEDED_BY_MANUAL_BUILD ${TARGETS_NEEDED_BY_MANUAL_BUILD} PARENT_SCOPE)
