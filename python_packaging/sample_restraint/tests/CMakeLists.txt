# Activate the CTest module for this directory.
enable_testing()


if(DOWNLOAD_GOOGLETEST)
    # Download and unpack googletest at configure time
    configure_file(CMakeLists.gtest.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
        message(WARNING "CMake step for googletest failed: ${result}")
        #    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
        message(WARNING "Build step for googletest failed: ${result}")
        #    message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                     EXCLUDE_FROM_ALL)
    add_library(GTest::Main ALIAS gtest_main)
else()
    find_package(GTest REQUIRED)
endif() # DOWNLOAD_GOOGLETEST

#
# Set up our tests
#


# Copy test files
# TODO: solve test data file requirement in a follow up change
file(DOWNLOAD
     https://github.com/kassonlab/sample_restraint/raw/master/tests/data/topol.tpr
     ${CMAKE_CURRENT_SOURCE_DIR}/topol.tpr
     STATUS _download_status
     LOG _download_log)
list(GET _download_status 0 _status)
if(${_status} GREATER 0)
    message(WARNING "Could not download test data: ${_download_log}")
else()
    configure_file(topol.tpr ${CMAKE_CURRENT_BINARY_DIR}/topol.tpr COPYONLY)
endif()
unset(_status)

# Note: Expects topol.tpr to be in CURRENT_BINARY_DUR
configure_file(testingconfiguration.in.h testingconfiguration.h)

# Test that the library can access its dependencies and build.
add_executable(library-test test_binding.cpp)
gtest_add_tests(TARGET library-test
                TEST_LIST BasicPlugin)

target_include_directories(library-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(library-test Gromacs::gmxapi gtest_main)

set_tests_properties(${BasicPlugin} PROPERTIES
                     COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}"
                     COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")


# Test the C++ code implementing the basic harmonic potential.
add_executable(harmonic-test test_harmonic.cpp)
gtest_add_tests(TARGET harmonic-test
                TEST_LIST HarmonicPotentialPlugin)

target_include_directories(harmonic-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(harmonic-test PROPERTIES SKIP_BUILD_RPATH FALSE)
target_link_libraries(harmonic-test harmonicpotential Gromacs::gmxapi
                      gtest_main)

set_tests_properties(${HarmonicPotentialPlugin} PROPERTIES
                     COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}"
                     COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")

# Test the C++ force evaluation for the restrained-ensemble biasing potential.
add_executable(histogram-test test_histogram.cpp)
gtest_add_tests(TARGET histogram-test
                TEST_LIST EnsembleHistogramPotentialPlugin)

target_include_directories(histogram-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(histogram-test PROPERTIES SKIP_BUILD_RPATH FALSE)
target_link_libraries(histogram-test ensemblepotential Gromacs::gmxapi
                      gtest_main)

set_tests_properties(${EnsembleHistogramPotentialPlugin} PROPERTIES
                     COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}"
                     COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")

# Test the flat-bottom bounding potential built in to the ensemble restraint.
add_executable(bounding-test test_bounding_restraint.cpp)
gtest_add_tests(TARGET bounding-test
                TEST_LIST EnsembleBoundingPotentialPlugin)

target_include_directories(bounding-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(bounding-test PROPERTIES SKIP_BUILD_RPATH FALSE)
target_link_libraries(bounding-test ensemblepotential Gromacs::gmxapi
                      GTest::Main)

set_tests_properties(${EnsembleBoundingPotentialPlugin} PROPERTIES
                     COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}"
                     COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
