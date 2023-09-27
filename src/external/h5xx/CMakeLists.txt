cmake_minimum_required(VERSION 3.1)

set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/platform.cmake")

# C language is only needed due to regression in FindHDF5.
project(H5XX LANGUAGES C CXX)

set(H5XX_VERSION 0.9.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(MPI "Build with MPI support for parallel IO" OFF)

# The FindBoost CMake module prefers multi-threaded libraries (filenames with
# postfix "-mt") over non-multi-threaded libraries. On Redhat or SuSE with
# installed system Boost libraries, this causes the system libraries (with
# "-mt") to override the custom-compiled libraries (without "-mt").
if(NOT DEFINED Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED FALSE CACHE BOOL "Prefer multi-threaded Boost libraries, with -mt suffix")
endif(NOT DEFINED Boost_USE_MULTITHREADED)

find_package(Boost 1.40.0 QUIET REQUIRED COMPONENTS unit_test_framework)
find_package(HDF5 COMPONENTS C REQUIRED)

# use
# $ cmake -DMPI=ON ...
# to configure with MPI support

if (MPI)
   find_package(MPI)
endif()
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
   link_libraries(${MPI_LIBRARIES})
   add_definitions(-DH5XX_USE_MPI -DMPICH_IGNORE_CXX_SEEK)
endif()

include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR})

#
# define target as a header-only library
#
add_library(h5xx INTERFACE)
target_compile_features(h5xx INTERFACE cxx_std_11)
target_include_directories(h5xx INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

include("${CMAKE_SOURCE_DIR}/cmake/install.cmake")

#
# provide unit tests
#
enable_testing()
include(CTest)
add_subdirectory(test)

