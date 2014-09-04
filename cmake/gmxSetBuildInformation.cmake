#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# - Check the username performing the build, as well as date and time
#
# gmx_set_build_information()
#
# The macro variables will be set to the user/host/cpu used for configuration,
# or anonymous/unknown if it cannot be detected (windows)
#
# BUILD_TIME
# BUILD_USER
# BUILD_HOST
# BUILD_CPU_VENDOR
# BUILD_CPU_BRAND
# BUILD_CPU_FAMILY
# BUILD_CPU_MODEL
# BUILD_CPU_STEPPING
# BUILD_CPU_FEATURES
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)

macro(gmx_set_build_information)
    IF(NOT DEFINED BUILD_USER)

    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM")
    else()
        set(GCC_INLINE_ASM_DEFINE "")
    endif()

    message(STATUS "Setting build user/date/host/cpu information")
    if(CMAKE_HOST_UNIX)
        execute_process( COMMAND date     OUTPUT_VARIABLE TMP_TIME    OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND whoami   OUTPUT_VARIABLE TMP_USER       OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND hostname OUTPUT_VARIABLE TMP_HOSTNAME   OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(BUILD_USER    "${TMP_USER}\@${TMP_HOSTNAME} [CMAKE]" CACHE INTERNAL "Build user")
        set(BUILD_TIME    "${TMP_TIME}" CACHE INTERNAL "Build date & time")
        execute_process( COMMAND uname -srm OUTPUT_VARIABLE TMP_HOST OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(BUILD_HOST    "${TMP_HOST}" CACHE INTERNAL "Build host & architecture")
        message(STATUS "Setting build user & time - OK")
    else()
        set(BUILD_USER    "Anonymous\@unknown [CMAKE]" CACHE INTERNAL "Build user")
        set(BUILD_TIME    "Unknown date" CACHE INTERNAL "Build date & time")
        set(BUILD_HOST    "${CMAKE_HOST_SYSTEM} ${CMAKE_HOST_SYSTEM_PROCESSOR}" CACHE INTERNAL "Build host & architecture")
        message(STATUS "Setting build user & time - not on Unix, using anonymous")
    endif()

    if(NOT CMAKE_CROSSCOMPILING)
        # Get CPU information, e.g. for deciding what SIMD support exists
        set(_compile_definitions "${GCC_INLINE_ASM_DEFINE} -I${CMAKE_SOURCE_DIR}/src -DGMX_CPUID_STANDALONE")
        if(GMX_TARGET_X86)
            set(_compile_definitions "${_compile_definitions} -DGMX_TARGET_X86")
        endif()
        try_run(GMX_CPUID_RUN_VENDOR GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_VENDOR ARGS "-vendor")
        try_run(GMX_CPUID_RUN_BRAND GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_BRAND ARGS "-brand")
        try_run(GMX_CPUID_RUN_FAMILY GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_FAMILY ARGS "-family")
        try_run(GMX_CPUID_RUN_MODEL GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_MODEL ARGS "-model")
       try_run(GMX_CPUID_RUN_STEPPING GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_STEPPING ARGS "-stepping")
        try_run(GMX_CPUID_RUN_FEATURES GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_FEATURES ARGS "-features")
        unset(_compile_definitions)

        string(STRIP "${OUTPUT_CPU_VENDOR}" OUTPUT_CPU_VENDOR)
        string(STRIP "${OUTPUT_CPU_BRAND}" OUTPUT_CPU_BRAND)
        string(STRIP "${OUTPUT_CPU_FAMILY}" OUTPUT_CPU_FAMILY)
        string(STRIP "${OUTPUT_CPU_MODEL}" OUTPUT_CPU_MODEL)
        string(STRIP "${OUTPUT_CPU_STEPPING}" OUTPUT_CPU_STEPPING)
        string(STRIP "${OUTPUT_CPU_FEATURES}" OUTPUT_CPU_FEATURES)

        if(GMX_CPUID_RUN_VENDOR EQUAL 0)
            set(BUILD_CPU_VENDOR   "${OUTPUT_CPU_VENDOR}"   CACHE INTERNAL "Build CPU vendor")
        else()
            set(BUILD_CPU_VENDOR   "Unknown, detect failed" CACHE INTERNAL "Build CPU vendor")
        endif()
        if(GMX_CPUID_RUN_BRAND EQUAL 0)
            set(BUILD_CPU_BRAND    "${OUTPUT_CPU_BRAND}"    CACHE INTERNAL "Build CPU brand")
        else()
            set(BUILD_CPU_BRAND    "Unknown, detect failed" CACHE INTERNAL "Build CPU brand")
        endif()
        if(GMX_CPUID_RUN_FAMILY EQUAL 0)
            set(BUILD_CPU_FAMILY   "${OUTPUT_CPU_FAMILY}"   CACHE INTERNAL "Build CPU family")
        else()
            set(BUILD_CPU_FAMILY   "0"                     CACHE INTERNAL "Build CPU family")
        endif()
        if(GMX_CPUID_RUN_MODEL EQUAL 0)
            set(BUILD_CPU_MODEL    "${OUTPUT_CPU_MODEL}"    CACHE INTERNAL "Build CPU model")
        else()
            set(BUILD_CPU_MODEL    "0"                     CACHE INTERNAL "Build CPU model")
        endif()
        if(GMX_CPUID_RUN_STEPPING EQUAL 0)
            set(BUILD_CPU_STEPPING "${OUTPUT_CPU_STEPPING}" CACHE INTERNAL "Build CPU stepping")
        else()
            set(BUILD_CPU_STEPPING "0"                     CACHE INTERNAL "Build CPU stepping")
        endif()
            if(GMX_CPUID_RUN_FEATURES EQUAL 0)
            set(BUILD_CPU_FEATURES "${OUTPUT_CPU_FEATURES}" CACHE INTERNAL "Build CPU features")
        else()
            set(BUILD_CPU_FEATURES ""                      CACHE INTERNAL "Build CPU features")
        endif()

    else()

        set(BUILD_CPU_VENDOR   "Unknown, cross-compiled"   CACHE INTERNAL "Build CPU vendor")
        set(BUILD_CPU_BRAND    "Unknown, cross-compiled"    CACHE INTERNAL "Build CPU brand")
        set(BUILD_CPU_FAMILY   "0"   CACHE INTERNAL "Build CPU family")
        set(BUILD_CPU_MODEL    "0"    CACHE INTERNAL "Build CPU model")
        set(BUILD_CPU_STEPPING "0" CACHE INTERNAL "Build CPU stepping")
        set(BUILD_CPU_FEATURES "" CACHE INTERNAL "Build CPU features")

    endif()

    ENDIF(NOT DEFINED BUILD_USER)
endmacro(gmx_set_build_information)
