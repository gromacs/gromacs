configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/config.h)

if (BUILD_TESTING)
    add_subdirectory(testutils)
    include(testutils/TestMacros.cmake)
endif (BUILD_TESTING)

add_subdirectory(gromacs)
add_subdirectory(kernel)
add_subdirectory(programs)

if(NOT GMX_FAHCORE)
  add_subdirectory(tools)
  add_subdirectory(ngmx)
endif(NOT GMX_FAHCORE)
