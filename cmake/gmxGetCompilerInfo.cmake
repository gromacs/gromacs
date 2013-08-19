#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#
# This macro attempts to parse the version string of the C compiler in use.
# With CMake 2.8.9 CMake provides a CMAKE_[C|CXX]_COMPILER_VERSION variable
# so we will use that if available.
#
# Currently supported are:
# - with cmake >2.8.8 all compilers supported by CMake
# - with cmake <=2.8.8: compilers that accept "-dumpversion" argument:
#   gcc, Intel Compiler (on Linux and Mac OS), Open64, EkoPath, clang
#   (and probably other gcc-compatible compilers).
# - with cmake <=2.8.8: xlC is not supported (it does not take -dumpversion,
#   but fortunately so far GROMACS never needs to know the version number)
#
# C_COMPILER_VERSION    - version string of the current C compiler (CMAKE_C_COMPILER)
# CXX_COMPILER_VERSION  - version string of the current C++ compiler (CMAKE_CXX_COMPILER)
#
macro(get_compiler_version)
    if(NOT C_COMPILER_VERSION)
        set(_cc_dumpversion_res 0)
        if (DEFINED CMAKE_C_COMPILER_VERSION AND CMAKE_VERSION VERSION_GREATER 2.8.8)
            set(_cc_version ${CMAKE_C_COMPILER_VERSION})
        elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
            set(_cc_dumpversion_res 1)
        else()
            execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                RESULT_VARIABLE _cc_dumpversion_res
                OUTPUT_VARIABLE _cc_version
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()

        if (${_cc_dumpversion_res} EQUAL 0)
            SET(C_COMPILER_VERSION ${_cc_version}
                CACHE STRING "C compiler version" FORCE)
        else ()
            SET(C_COMPILER_VERSION ""
                CACHE STRING "C compiler version not available" FORCE)
        endif ()
    endif()

    if(NOT CXX_COMPILER_VERSION AND CMAKE_CXX_COMPILER_LOADED)
        set(_cxx_dumpversion_res 0)
        if (DEFINED CMAKE_CXX_COMPILER_VERSION AND CMAKE_VERSION VERSION_GREATER 2.8.8)
            set(_cxx_version ${CMAKE_CXX_COMPILER_VERSION})
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
            set(_cxx_dumpversion_res 1)
        else()
            execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                RESULT_VARIABLE _cxx_dumpversion_res
                OUTPUT_VARIABLE _cxx_version
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()

        if (${_cxx_dumpversion_res} EQUAL 0)
            SET(CXX_COMPILER_VERSION ${_cxx_version}
                CACHE STRING "C++ compiler version string" FORCE)
        else ()
            SET(CXX_COMPILER_VERSION ""
                CACHE STRING "C++ compiler version string not available" FORCE)
        endif ()
    endif ()

    if (NOT "${C_COMPILER_VERSION}" STREQUAL "${CXX_COMPILER_VERSION}" AND CMAKE_CXX_COMPILER_LOADED)
        message(WARNING "The version of the C and C++ compilers does not match. Note that mixing different C/C++ compilers can cause problems!")
    endif ()

    mark_as_advanced(C_COMPILER_VERSION CXX_COMPILER_VERSION)
endmacro()

# This macro attempts to get a reasonable version string for a compiler,
# and also extracts compiler flags.
#
# Parameters:
#   LANGUAGE       - C or CXX, the compiler to check for
#   BUILD_COMPILER - [output variable] string with compiler path, ID and
#                    some compiler-provided information
#   BUILD_FLAGS    - [output variable] flags for the compiler
#
macro(get_compiler_info LANGUAGE BUILD_COMPILER BUILD_FLAGS)
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        set(_flag_to_query_version "-qversion")
    else()
        set(_flag_to_query_version "--version")
    endif()
    execute_process(COMMAND ${CMAKE_${LANGUAGE}_COMPILER} ${_flag_to_query_version}
        RESULT_VARIABLE _exec_result
        OUTPUT_VARIABLE _compiler_version
        ERROR_VARIABLE  _compiler_version)

    if(_exec_result)
        # Try executing just the compiler command, since --version failed
        execute_process(COMMAND ${CMAKE_${LANGUAGE}_COMPILER}
            RESULT_VARIABLE _exec_result
            OUTPUT_VARIABLE _compiler_version
            ERROR_VARIABLE  _compiler_version)
    endif()
    if(NOT "${_compiler_version}" STREQUAL "")
        string(REGEX MATCH "[^\n]+" _compiler_version "${_compiler_version}")
    endif()

    set(${BUILD_COMPILER}
        "${CMAKE_${LANGUAGE}_COMPILER} ${CMAKE_${LANGUAGE}_COMPILER_ID} ${_compiler_version}")
    set(_build_flags "${CMAKE_${LANGUAGE}_FLAGS}")
    string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
    set(_build_flags "${_build_flags} ${CMAKE_${LANGUAGE}_FLAGS_${_build_type}}")
    set(${BUILD_FLAGS} ${_build_flags})
endmacro()
