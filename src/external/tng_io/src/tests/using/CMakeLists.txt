# Test project to test building with libtng_io
#
# Combinations to test for full coverage:
#   TEST_ZLIB=ON (22 combinations)
#     ZLIB_LIBRARY=shared/static
#     TEST_BUNDLING=ON (4 combinations)
#       BUILD_SHARED_LIBS=ON/OFF
#       TEST_EXECUTABLE=ON/OFF
#     TEST_BUNDLING=OFF (7 combinations)
#       TEST_EXECUTABLE=ON
#         TNG_IO_DIR=sh.sh/st.sh/st.st
#       TEST_EXECUTABLE=OFF
#         BUILD_SHARED_LIBS=ON
#           TNG_IO_DIR=sh.sh/st.sh/st.st
#         BUILD_SHARED_LIBS=OFF
#           TNG_IO_DIR=sh.sh
#   TEST_ZLIB=OFF (27 combinations)
#     TEST_BUNDLING=ON (12 combinations)
#       ZLIB_LIBRARY=shared/static/none
#       BUILD_SHARED_LIBS=ON/OFF
#       TEST_EXECUTABLE=ON/OFF
#     TEST_BUNDLING=OFF (15 combinations)
#       TEST_EXECUTABLE=ON
#         TNG_IO_DIR=sh.sh/sh.st/sh.int/st.sh/st.st/st.int
#       TEST_EXECUTABLE=OFF
#         BUILD_SHARED_LIBS=ON
#           TNG_IO_DIR=sh.sh/sh.st/sh.int/st.sh/st.st/st.int
#         BUILD_SHARED_LIBS=OFF
#           TNG_IO_DIR=sh.sh/sh.st/sh.int
#
# Combinations that cannot work:
#   TEST_ZLIB=ON, tng built with internal zlib
#   TEST_ZLIB=ON, tng built with BUILD_SHARED_LIBS=ON + static zlib
#   BUILD_SHARED_LIBS=ON, tng built as static (unless compiled with PIC)

cmake_minimum_required(VERSION 3.1)

project(tng_io_test)

option(BUILD_SHARED_LIBS "Test building a shared library" ON)
option(TEST_EXECUTABLE "Test building an executable instead of a library" OFF)
option(TEST_BUNDLING "Test bundling tng" OFF)
option(TEST_ZLIB "Test with zlib in using code" OFF)
option(TEST_BUNDLED_ZLIB "Test bundling zlib with tng" OFF)
if (NOT TEST_BUNDLING)
    set(TEST_BUNDLED_ZLIB OFF)
endif()

if (TEST_ZLIB OR NOT TEST_BUNDLED_ZLIB)
    find_package(ZLIB REQUIRED)
endif()

if (TEST_BUNDLING)
    include(../../../BuildTNG.cmake)
    if (TEST_BUNDLED_ZLIB)
        message(STATUS "Bundling tng_io with internal ZLIB")
        add_tng_io_library(tng_io OBJECT OWN_ZLIB)
    else()
        message(STATUS "Bundling tng_io with external ZLIB")
        add_tng_io_library(tng_io OBJECT)
    endif()
    add_library(tng_io::tng_io ALIAS tng_io)
else()
    message(STATUS "Using external tng_io")
    find_package(TNG_IO REQUIRED)
endif()

set(SOURCES use_tng.c)
if (TEST_ZLIB)
    message(STATUS "Using zlib in test application")
    list(APPEND SOURCES use_zlib.c)
else()
    message(STATUS "Not using zlib in test application")
    list(APPEND SOURCES dummy_zlib.c)
endif()

if (TEST_EXECUTABLE)
    message(STATUS "Building an executable linked against tng_io")
    add_executable(test_target main.c ${SOURCES})
else()
    message(STATUS "Building a library linked against tng_io")
    add_library(test_target ${SOURCES})
    install(TARGETS test_target EXPORT test DESTINATION lib)
    install(EXPORT test DESTINATION lib/cmake/test)
    add_executable(test_exe main.c)
    target_link_libraries(test_exe PRIVATE test_target)
endif()

if (TEST_BUNDLING)
    target_link_libraries(test_target PRIVATE $<BUILD_INTERFACE:tng_io::tng_io>)
else()
    target_link_libraries(test_target PRIVATE tng_io::tng_io)
endif()
if (TEST_ZLIB)
    target_link_libraries(test_target PRIVATE ZLIB::ZLIB)
endif()
