configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/config.h)

if (BUILD_TESTING)
    include(testutils/TestMacros.cmake)
    add_subdirectory(testutils)
endif (BUILD_TESTING)

add_subdirectory(gromacs)
add_subdirectory(kernel)
add_subdirectory(programs)

if(NOT GMX_FAHCORE)
  add_subdirectory(tools)
  add_subdirectory(ngmx)
  add_subdirectory(contrib)
endif(NOT GMX_FAHCORE)

if(CPPCHECK_FOUND)
  include(../cmake/CppcheckTargets.cmake)
  file(GLOB CPPCHECK_SOURCE_FILES */*.cpp */*/*.cpp */*/*/*.cpp)
  add_cppcheck_sources(cppcheck UNUSED_FUNCTIONS STYLE ${CPPCHECK_SOURCE_FILES} -DHAVE_CONFIG_H -DCPPCHECK -Igromacs/legacyheaders/)
endif()
