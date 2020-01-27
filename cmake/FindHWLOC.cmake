#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
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

###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find HWLOC include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HWLOC
#               [REQUIRED]) # Fail with error if hwloc is not found
#
# This module finds headers and hwloc library.
# Results are reported in variables:
#  HWLOC_FOUND           - True if headers and requested libraries were found
#  HWLOC_INCLUDE_DIRS    - hwloc include directories
#  HWLOC_LIBRARY_DIRS    - Link directories for hwloc libraries
#  HWLOC_LIBRARIES       - hwloc component libraries to be linked
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DHWLOC_DIR=path/to/hwloc):
#  HWLOC_DIR             - Where to find the base directory of hwloc
#  HWLOC_INCDIR          - Where to find the header files
#  HWLOC_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: HWLOC_DIR, HWLOC_INCDIR, HWLOC_LIBDIR

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)

include(CheckStructHasMember)
include(CheckCSourceCompiles)

if (NOT HWLOC_FOUND)
  set(HWLOC_DIR "" CACHE PATH "Installation directory of HWLOC library")
endif()

set(ENV_HWLOC_DIR "$ENV{HWLOC_DIR}")
set(ENV_HWLOC_INCDIR "$ENV{HWLOC_INCDIR}")
set(ENV_HWLOC_LIBDIR "$ENV{HWLOC_LIBDIR}")
set(HWLOC_GIVEN_BY_USER "FALSE")
if ( HWLOC_DIR OR ( HWLOC_INCDIR AND HWLOC_LIBDIR) OR ENV_HWLOC_DIR OR (ENV_HWLOC_INCDIR AND ENV_HWLOC_LIBDIR) )
  set(HWLOC_GIVEN_BY_USER "TRUE")
endif()

if (NOT HWLOC_FIND_QUIETLY)
    message(STATUS "Looking for HWLOC")
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
if(ENV_HWLOC_INCDIR)
    list(APPEND _inc_env "${ENV_HWLOC_INCDIR}")
elseif(ENV_HWLOC_DIR)
    list(APPEND _inc_env "${ENV_HWLOC_DIR}")
    list(APPEND _inc_env "${ENV_HWLOC_DIR}/include")
    list(APPEND _inc_env "${ENV_HWLOC_DIR}/include/hwloc")
else()
    if(WIN32)
        string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
    else()
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
    endif()
endif()
list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)

# set paths where to look for
set(PATH_TO_LOOK_FOR "${_inc_env}")

# Try to find the hwloc header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(HWLOC_INCDIR)
    set(HWLOC_hwloc.h_DIRS "HWLOC_hwloc.h_DIRS-NOTFOUND")
    find_path(HWLOC_hwloc.h_DIRS
        NAMES hwloc.h
        HINTS ${HWLOC_INCDIR})
else()
    if(HWLOC_DIR)
        set(HWLOC_hwloc.h_DIRS "HWLOC_hwloc.h_DIRS-NOTFOUND")
        find_path(HWLOC_hwloc.h_DIRS
	    NAMES hwloc.h
	    HINTS ${HWLOC_DIR}
	    PATH_SUFFIXES "include" "include/hwloc")
    else()
        set(HWLOC_hwloc.h_DIRS "HWLOC_hwloc.h_DIRS-NOTFOUND")
        find_path(HWLOC_hwloc.h_DIRS
	    NAMES hwloc.h
	    HINTS ${PATH_TO_LOOK_FOR}
	    PATH_SUFFIXES "hwloc")
    endif()
endif()
mark_as_advanced(HWLOC_hwloc.h_DIRS)

# Add path to cmake variable
# ------------------------------------
if (HWLOC_hwloc.h_DIRS)
    set(HWLOC_INCLUDE_DIRS "${HWLOC_hwloc.h_DIRS}")
else ()
    set(HWLOC_INCLUDE_DIRS "HWLOC_INCLUDE_DIRS-NOTFOUND")
    if(NOT HWLOC_FIND_QUIETLY)
        message(STATUS "Looking for hwloc -- hwloc.h not found")
    endif()
endif ()

if (HWLOC_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES HWLOC_INCLUDE_DIRS)
endif ()


# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
if(ENV_HWLOC_LIBDIR)
    list(APPEND _lib_env "${ENV_HWLOC_LIBDIR}")
elseif(ENV_HWLOC_DIR)
    list(APPEND _lib_env "${ENV_HWLOC_DIR}")
    list(APPEND _lib_env "${ENV_HWLOC_DIR}/lib")
else()
    if(WIN32)
        string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
    else()
        if(APPLE)
	    string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
        else()
	    string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
        endif()
        list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
        list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    endif()
endif()
list(REMOVE_DUPLICATES _lib_env)

# set paths where to look for
set(PATH_TO_LOOK_FOR "${_lib_env}")

# Try to find the hwloc lib in the given paths
# ----------------------------------------------

# call cmake macro to find the lib path
if(HWLOC_LIBDIR)
    set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")
    find_library(HWLOC_hwloc_LIBRARY
        NAMES hwloc
        HINTS ${HWLOC_LIBDIR})
else()
    if(HWLOC_DIR)
        set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")
        find_library(HWLOC_hwloc_LIBRARY
	    NAMES hwloc
	    HINTS ${HWLOC_DIR}
	    PATH_SUFFIXES lib lib32 lib64)
    else()
        set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")
        find_library(HWLOC_hwloc_LIBRARY
	    NAMES hwloc
	    HINTS ${PATH_TO_LOOK_FOR})
    endif()
endif()
mark_as_advanced(HWLOC_hwloc_LIBRARY)

# If found, add path to cmake variable
# ------------------------------------
if (HWLOC_hwloc_LIBRARY)
    get_filename_component(hwloc_lib_path ${HWLOC_hwloc_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(HWLOC_LIBRARIES    "${HWLOC_hwloc_LIBRARY}")
    set(HWLOC_LIBRARY_DIRS "${hwloc_lib_path}")
else ()
    set(HWLOC_LIBRARIES    "HWLOC_LIBRARIES-NOTFOUND")
    set(HWLOC_LIBRARY_DIRS "HWLOC_LIBRARY_DIRS-NOTFOUND")
    if(NOT HWLOC_FIND_QUIETLY)
        message(STATUS "Looking for hwloc -- lib hwloc not found")
    endif()
endif ()

if (HWLOC_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES HWLOC_LIBRARY_DIRS)
endif ()

# check a function to validate the find
if(HWLOC_LIBRARIES)

    set(REQUIRED_INCDIRS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # HWLOC
    if (HWLOC_INCLUDE_DIRS)
        set(REQUIRED_INCDIRS "${HWLOC_INCLUDE_DIRS}")
    endif()
    if (HWLOC_LIBRARY_DIRS)
        set(REQUIRED_LIBDIRS "${HWLOC_LIBRARY_DIRS}")
    endif()
    set(REQUIRED_LIBS "${HWLOC_LIBRARIES}")

    # set required libraries for link
    set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
    set(CMAKE_REQUIRED_LIBRARIES)
    foreach(lib_dir ${REQUIRED_LIBDIRS})
        list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
    endforeach()
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
    string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

    # test link
    include(CheckFunctionExists)
    check_function_exists(hwloc_topology_init HWLOC_WORKS)
    mark_as_advanced(HWLOC_WORKS)

    if(NOT HWLOC_WORKS)
        if(NOT HWLOC_FIND_QUIETLY)
	    message(STATUS "Looking for hwloc : test of hwloc_topology_init with hwloc library fails")
	    message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
	    message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
	    message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
        endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
endif(HWLOC_LIBRARIES)

if (HWLOC_LIBRARIES)
  if (HWLOC_LIBRARY_DIRS)
    list(GET HWLOC_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET HWLOC_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(HWLOC_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of HWLOC library" FORCE)
  else()
    set(HWLOC_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of HWLOC library" FORCE)
  endif()
endif()
mark_as_advanced(HWLOC_DIR)
mark_as_advanced(HWLOC_DIR_FOUND)

# Function for converting hex version numbers from HWLOC_API_VERSION if necessary
function(HEX2DEC str res)
    string(LENGTH "${str}" len)
    if("${len}" EQUAL 1)
        if("${str}" MATCHES "[0-9]")
            set(${res} "${str}" PARENT_SCOPE)
        elseif( "${str}" MATCHES "[aA]")
            set(${res} 10 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[bB]")
            set(${res} 11 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[cC]")
            set(${res} 12 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[dD]")
            set(${res} 13 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[eE]")
            set(${res} 14 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[fF]")
            set(${res} 15 PARENT_SCOPE)
        else()
            return()
        endif()
    else()
        string(SUBSTRING "${str}" 0 1 str1)
        string(SUBSTRING "${str}" 1 -1 str2)
        hex2dec(${str1} res1)
        hex2dec(${str2} res2)
        math(EXPR val "16 * ${res1} + ${res2}")
        set(${res} "${val}" PARENT_SCOPE)
    endif()
endfunction()

if(HWLOC_INCLUDE_DIRS)
    # If we are not cross-compiling we try to use the hwloc-info program
    if(NOT CMAKE_CROSSCOMPILING)
        find_program(HWLOC_INFO "hwloc-info"
            HINTS "${HWLOC_DIR}" ENV HWLOC_DIR
            PATH_SUFFIXES bin
            )
        mark_as_advanced(HWLOC_INFO)

        if(HWLOC_INFO)
            execute_process(COMMAND ${HWLOC_INFO} "--version"
                            RESULT_VARIABLE HWLOC_INFO_RES
                            OUTPUT_VARIABLE HWLOC_INFO_OUT
                            ERROR_VARIABLE  HWLOC_INFO_ERR)

            if(HWLOC_INFO_ERR)
	        message(STATUS "Error executing hwloc-info: ${HWLOC_INFO_ERR}")
            endif()
            string(REGEX MATCH "[0-9]+.*[0-9]+" HWLOC_INFO_OUT "${HWLOC_INFO_OUT}")
            set(HWLOC_LIBRARY_VERSION ${HWLOC_INFO_OUT})
        endif()
    endif()

    if (NOT HWLOC_FIND_QUIETLY)
        message(STATUS "hwloc version: ${HWLOC_VERSION}")
    endif()

    # Parse header if cross-compiling, or if hwloc-info was not found
    # Also used to check that library and header versions match
    # HWLOC is never installed as a framework on OS X, so this should always work.
    file(READ "${HWLOC_INCLUDE_DIRS}/hwloc.h"
         HEADER_CONTENTS LIMIT 16384)
    string(REGEX REPLACE ".*#define HWLOC_API_VERSION (0[xX][0-9a-fA-F]+).*" "\\1"
           HWLOC_API_VERSION "${HEADER_CONTENTS}")
    string(SUBSTRING "${HWLOC_API_VERSION}" 4 2 HEX_MAJOR)
    string(SUBSTRING "${HWLOC_API_VERSION}" 6 2 HEX_MINOR)
    string(SUBSTRING "${HWLOC_API_VERSION}" 8 2 HEX_PATCH)
    hex2dec(${HEX_MAJOR} DEC_MAJOR)
    hex2dec(${HEX_MINOR} DEC_MINOR)
    hex2dec(${HEX_PATCH} DEC_PATCH)
    set(HWLOC_HEADER_VERSION "${DEC_MAJOR}.${DEC_MINOR}.${DEC_PATCH}")

    if (HWLOC_LIBRARY_VERSION AND HWLOC_HEADER_VERSION)
        string(SUBSTRING "${HWLOC_LIBRARY_VERSION}" 0 1 LIBRARY_MAJOR)
        string(SUBSTRING "${HWLOC_HEADER_VERSION}" 0 1 HEADER_MAJOR)
        string(COMPARE EQUAL "${LIBRARY_MAJOR}" "${HEADER_MAJOR}" HWLOC_VERSION_CHECK)
        if(NOT HWLOC_VERSION_CHECK)
            message(FATAL_ERROR "Detected version mismatch between HWLOC headers and library. "
            "Library version is ${HWLOC_LIBRARY_VERSION}, but header version is ${HWLOC_HEADER_VERSION}. "
            "Make sure that you have the correct include and library directory set for HWLOC")
        endif()
    endif()
    if (HWLOC_LIBRARY_VERSION)
        set(HWLOC_VERSION ${HWLOC_LIBRARY_VERSION} CACHE STRING "HWLOC library version")

    else()
        set(HWLOC_VERSION ${HWLOC_HEADER_VERSION} CACHE STRING "HWLOC library version")
    endif()
    set(GMX_HWLOC_API_VERSION ${HWLOC_API_VERSION} CACHE STRING "HWLOC API version during configuration time")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC
                                  REQUIRED_VARS HWLOC_LIBRARIES HWLOC_INCLUDE_DIRS
                                  VERSION_VAR HWLOC_VERSION)

mark_as_advanced(HWLOC_INCLUDE_DIRS HWLOC_LIBRARIES HWLOC_VERSION)
