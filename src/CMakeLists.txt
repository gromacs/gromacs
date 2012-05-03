configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/config.h)

if (BUILD_TESTING)
    include(testutils/TestMacros.cmake)
    add_subdirectory(testutils)
endif (BUILD_TESTING)

add_subdirectory(gromacs)
add_subdirectory(programs)

if(NOT GMX_FAHCORE)
  add_subdirectory(tools)
  add_subdirectory(ngmx)
  add_subdirectory(contrib)
endif(NOT GMX_FAHCORE)
