#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2009- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    message(FATAL_ERROR "CMake called with the wrong directory; please use the top GROMACS directory (${CMAKE_CURRENT_SOURCE_DIR}/../)")
endif()

######################################
# Output compiler and CFLAGS used
######################################
include(GetCompilerInfo.cmake)
get_compiler_info(C BUILD_C_COMPILER)
get_compiler_info(CXX BUILD_CXX_COMPILER)
if(GMX_GPU_CUDA)
    get_cuda_compiler_info(CUDA_COMPILER_INFO CUDA_DEVICE_COMPILER_FLAGS CUDA_HOST_COMPILER_FLAGS)
endif()

if(GMX_GPU_HIP)
    get_hip_compiler_info(HIP_COMPILER_INFO HIP_COMPILER_FLAGS)
endif()

# Make a file with compiler flags used for libgromacs for each
# langauge and build configuration.  The one that corresponds to
# CMAKE_BUILD_TYPE is #included into buildinfo.h and populates the
# fields e.g. printed to the log file.
file(GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/include/compilerflaginfo-$<CONFIG>-$<COMPILE_LANGUAGE>.h
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/compilerflaginfo.h.cmakein
    CONDITION $<CONFIG:${CMAKE_BUILD_TYPE}>
    )

####

option(GMX_CLANG_TIDY "Use clang-tidy" OFF)
mark_as_advanced(GMX_CLANG_TIDY)
if (GMX_CLANG_TIDY)
   if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
   elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithAssert")
   elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
   elseif("${CMAKE_BUILD_TYPE}" STREQUAL "ASAN")
   else()
       message(FATAL_ERROR "Can only use clang-tidy with build type containing asserts: Debug, RelWithAssert, RelWithDebInfo, ASAN.")
   endif()
   set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
   mark_as_advanced(CMAKE_EXPORT_COMPILE_COMMANDS)
   set(CLANG_TIDY "clang-tidy" CACHE STRING "Name of clang-tidy executable")
   find_program(CLANG_TIDY_EXE NAMES "${CLANG_TIDY}"
       DOC "Path to clang-tidy executable")
   if(NOT CLANG_TIDY_EXE)
       message(FATAL_ERROR "clang-tidy not found.")
   endif()
   mark_as_advanced(CLANG_TIDY)
   mark_as_advanced(CLANG_TIDY_EXE)
endif()

# Create a basic target for the `src` section of the build tree to capture
# the library-level shared details through CMake infrastructure. It is not
# installed or exported, so it must only be used as a PRIVATE dependency by
# installed targets.
# Initially, this is just an INTERFACE target to provide include directory.
# It should also absorb global variables and compiler/linker details to be
# provided as transitive usage requirements.
# It could expand to aggregate the module targets in the future.
add_library(common INTERFACE)
target_include_directories(common INTERFACE
                           ${CMAKE_CURRENT_SOURCE_DIR}/include
                           ${CMAKE_CURRENT_BINARY_DIR}/include)

add_subdirectory(external)

if (BUILD_TESTING)
    if(NOT GMX_DEVELOPER_BUILD)
        set(UNITTEST_TARGET_OPTIONS EXCLUDE_FROM_ALL)
    endif()
    include(testutils/TestMacros.cmake)
    add_subdirectory(testutils)
endif()

add_subdirectory(gromacs)
add_subdirectory(programs)

# Configure header files with configuration-specific values. This step
# should follow all introspection e.g. looking for headers and
# libraries. If not, cmake will need to change the contents of the
# file upon subsequent runs of cmake. This can mean that
#
#  cmake $src && make && make test
#
# requires building all the source files that depend on the changed
# header file in both of the make stages. That's slow, and is useless
# busy work for ccache, too.
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
configure_file(config.h.cmakein include/config.h)
configure_file(gmxpre-config.h.cmakein include/gmxpre-config.h)

set(CMAKE_BUILD_CONFIGURATION_C_FLAGS   ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
set(CMAKE_BUILD_CONFIGURATION_CXX_FLAGS ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
configure_file(buildinfo.h.cmakein include/buildinfo.h ESCAPE_QUOTES)
configure_file(contributors.h.cmakein include/contributors.h)
