#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013 by the GROMACS development team.
# Copyright (c) 2014,2015,2016,2017,2018 by the GROMACS development team.
# Copyright (c) 2019,2020, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

######################################
# Output compiler and CFLAGS used
######################################
include(GetCompilerInfo.cmake)
get_compiler_info(C BUILD_C_COMPILER)
get_compiler_info(CXX BUILD_CXX_COMPILER)
if(GMX_USE_CUDA)
    if(NOT GMX_CLANG_CUDA)
        GMX_SET_CUDA_NVCC_FLAGS()
    endif()

    get_cuda_compiler_info(CUDA_COMPILER_INFO CUDA_DEVICE_COMPILER_FLAGS CUDA_HOST_COMPILER_FLAGS)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
configure_file(config.h.cmakein config.h)
configure_file(gmxpre-config.h.cmakein gmxpre-config.h)
configure_file(buildinfo.h.cmakein buildinfo.h ESCAPE_QUOTES)
# Make a file with compiler flags used for libgromacs for each
# langauge and build configuration.  The one that corresponds to
# CMAKE_BUILD_TYPE is #included into buildinfo.h and populates the
# fields e.g. printed to the log file.
file(GENERATE
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/compilerflaginfo-$<CONFIG>-$<COMPILE_LANGUAGE>.h
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/compilerflaginfo.h.cmakein
    CONDITION $<CONFIG:${CMAKE_BUILD_TYPE}>
    )

####
set(IGNORED_CLANG_ALL_WARNINGS
    "-Wno-c++98-compat -Wno-c++98-compat-pedantic" #No intention of C++98 compability
    "-Wno-source-uses-openmp" #Don't warn for no-omp build
    "-Wno-c++17-extensions"   #Allowed in attributes (compilers are required to ignore unknown attributes)
    "-Wno-documentation-unknown-command" #Custom commands are used
    "-Wno-covered-switch-default" #GCC gives maybe-uninitialized without default label and checks for illegal enum values.
    "-Wno-switch-enum" # default statement for enum is OK

    # We need to use macros like
    # GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR. Those will look strange
    # if they don't have a semicolon after them, and might confuse
    # tools like IDEs also.
    "-Wno-extra-semi-stmt"

    #Following ones are undecided/TODO
    "-Wno-disabled-macro-expansion"
    "-Wno-cast-align"
    "-Wno-reserved-id-macro"
    "-Wno-global-constructors"
    "-Wno-exit-time-destructors"
    "-Wno-unused-macros"
    "-Wno-weak-vtables"
    "-Wno-conditional-uninitialized"
    "-Wno-format-nonliteral"
    "-Wno-shadow"
    "-Wno-cast-qual"
    "-Wno-documentation"
    "-Wno-used-but-marked-unused"
    "-Wno-padded"
    "-Wno-float-equal"
    "-Wno-old-style-cast"
    "-Wno-conversion"
    "-Wno-double-promotion")
string(REPLACE " " ";" IGNORED_CLANG_ALL_WARNINGS "${IGNORED_CLANG_ALL_WARNINGS}")

option(GMX_CLANG_TIDY "Use clang-tidy" OFF)
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

# Activate targets for new C++ API components and docs.
if(GMX_NATIVE_WINDOWS OR GMX_BUILD_MDRUN_ONLY OR NOT BUILD_SHARED_LIBS)
    # GMXAPI has not been tested in Microsoft environments.
    # GMXAPI relies on libgromacs and is incompatible with an `mdrun`-only build.
    # GMXAPI requires position-independent code
    set(_GMXAPI_DEFAULT OFF)
else()
    set(_GMXAPI_DEFAULT ON)
endif()
option(GMXAPI "Experimental external interfaces" ${_GMXAPI_DEFAULT})
if (GMXAPI)
   if (GMX_BUILD_MDRUN_ONLY)
       message(FATAL_ERROR "GMXAPI relies on libgromacs and is incompatible with GMX_BUILD_MDRUN_ONLY.")
   endif()
   if(NOT ${BUILD_SHARED_LIBS})
       # Note: this conditional should check for the existence of a libgromacs target supporting PIC
       # using the POSITION_INDEPENDENT_CODE property, but for now the only facility we have is the global
       # variable, BUILD_SHARED_LIBS.
       # TODO: gmxapi should gracefully build for dynamic linking or static linking for PIC or without.
       message(FATAL_ERROR "GMXAPI requires position-independent code. Set -DGMXAPI=OFF or -DBUILD_SHARED_LIBS=ON.")
   endif()
   add_subdirectory(api)
endif()
