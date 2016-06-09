#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM=1")
    else()
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM=0")
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

    # Set up some defaults that will usually be overridden
    if(CMAKE_CROSSCOMPILING)
        set(_reason ", cross-compiled")
    endif()
    set(OUTPUT_CPU_VENDOR   "Unknown${_reason}")
    set(OUTPUT_CPU_BRAND    "Unknown${_reason}")
    set(OUTPUT_CPU_FAMILY   "0")
    set(OUTPUT_CPU_MODEL    "0")
    set(OUTPUT_CPU_STEPPING "0")
    set(OUTPUT_CPU_FEATURES "Unknown${_reason}")
    unset(_reason)

    if(NOT CMAKE_CROSSCOMPILING)
        # Get CPU information, e.g. for deciding what SIMD support probably exists
        set(_compile_definitions "${GCC_INLINE_ASM_DEFINE} -I${CMAKE_SOURCE_DIR}/src -DGMX_CPUINFO_STANDALONE ${GMX_STDLIB_CXX_FLAGS}")
        if(GMX_TARGET_X86)
            set(_compile_definitions "${_compile_definitions} -DGMX_TARGET_X86")
        endif()

        set(GMX_BUILDINFORMATION_BINARY "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/GmxBuildInformation${CMAKE_EXECUTABLE_SUFFIX}")
        set(LINK_LIBRARIES "${GMX_STDLIB_LIBRARIES}")
        # TODO Extract this try_compile to a helper function, because
        # it duplicates code in gmxDetectSimd.cmake
        try_compile(GMX_BUILDINFORMATION_COMPILED
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/gromacs/hardware/cpuinfo.cpp"
            COMPILE_DEFINITIONS "${_compile_definitions}"
            CMAKE_FLAGS "-DLINK_LIBRARIES=${LINK_LIBRARIES}"
            OUTPUT_VARIABLE GMX_BUILDINFORMATION_COMPILED_OUTPUT
            COPY_FILE ${GMX_BUILDINFORMATION_BINARY})
        unset(_compile_definitions)

        if(GMX_BUILDINFORMATION_COMPILED)
            # TODO Extract this duplication to a helper function (also
            # from gmxDetectSimd.cmake)
            if(NOT DEFINED GMX_BUILDINFORMATION_RUN_VENDOR)
                execute_process(COMMAND ${GMX_BUILDINFORMATION_BINARY} "-vendor"
                    RESULT_VARIABLE GMX_BUILDINFORMATION_RUN_VENDOR
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_BUILDINFORMATION_RUN_VENDOR "${GMX_BUILDINFORMATION_RUN_VENDOR}" CACHE INTERNAL "Result of running cpuinfo code with arg -vendor")
                if(GMX_BUILDINFORMATION_RUN_VENDOR EQUAL 0)
                    string(STRIP "${OUTPUT_TMP}" OUTPUT_CPU_VENDOR)
                endif()
            endif()
            if(NOT DEFINED GMX_BUILDINFORMATION_RUN_BRAND)
                execute_process(COMMAND ${GMX_BUILDINFORMATION_BINARY} "-brand"
                    RESULT_VARIABLE GMX_BUILDINFORMATION_RUN_BRAND
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_BUILDINFORMATION_RUN_BRAND "${GMX_BUILDINFORMATION_RUN_BRAND}" CACHE INTERNAL "Result of running cpuinfo code with arg -brand")
                if(GMX_BUILDINFORMATION_RUN_BRAND EQUAL 0)
                    string(STRIP "${OUTPUT_TMP}" OUTPUT_CPU_BRAND)
                endif()
            endif()
            if(NOT DEFINED GMX_BUILDINFORMATION_RUN_FAMILY)
                execute_process(COMMAND ${GMX_BUILDINFORMATION_BINARY} "-family"
                    RESULT_VARIABLE GMX_BUILDINFORMATION_RUN_FAMILY
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_BUILDINFORMATION_RUN_FAMILY "${GMX_BUILDINFORMATION_RUN_FAMILY}" CACHE INTERNAL "Result of running cpuinfo code with arg -family")
                if(GMX_BUILDINFORMATION_RUN_FAMILY EQUAL 0)
                    string(STRIP "${OUTPUT_TMP}" OUTPUT_CPU_FAMILY)
                endif()
            endif()
            if(NOT DEFINED GMX_BUILDINFORMATION_RUN_MODEL)
                execute_process(COMMAND ${GMX_BUILDINFORMATION_BINARY} "-model"
                    RESULT_VARIABLE GMX_BUILDINFORMATION_RUN_MODEL
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_BUILDINFORMATION_RUN_MODEL "${GMX_BUILDINFORMATION_RUN_MODEL}" CACHE INTERNAL "Result of running cpuinfo code with arg -model")
                if(GMX_BUILDINFORMATION_RUN_MODEL EQUAL 0)
                    string(STRIP "${OUTPUT_TMP}" OUTPUT_CPU_MODEL)
                endif()
            endif()
            if(NOT DEFINED GMX_BUILDINFORMATION_RUN_STEPPING)
                execute_process(COMMAND ${GMX_BUILDINFORMATION_BINARY} "-stepping"
                    RESULT_VARIABLE GMX_BUILDINFORMATION_RUN_STEPPING
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_BUILDINFORMATION_RUN_STEPPING "${GMX_BUILDINFORMATION_RUN_STEPPING}" CACHE INTERNAL "Result of running cpuinfo code with arg -stepping")
                if(GMX_BUILDINFORMATION_RUN_STEPPING EQUAL 0)
                    string(STRIP "${OUTPUT_TMP}" OUTPUT_CPU_STEPPING)
                endif()
            endif()
            if(NOT DEFINED GMX_BUILDINFORMATION_RUN_FEATURES)
                execute_process(COMMAND ${GMX_BUILDINFORMATION_BINARY} "-features"
                    RESULT_VARIABLE GMX_BUILDINFORMATION_RUN_FEATURES
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_BUILDINFORMATION_RUN_FEATURES "${GMX_BUILDINFORMATION_RUN_FEATURES}" CACHE INTERNAL "Result of running cpuinfo code with arg -features")
                if(GMX_BUILDINFORMATION_RUN_FEATURES EQUAL 0)
                    string(STRIP "${OUTPUT_TMP}" OUTPUT_CPU_FEATURES)
                endif()
            endif()
        endif()
    endif()

    set(BUILD_CPU_VENDOR   "${OUTPUT_CPU_VENDOR}"   CACHE INTERNAL "Build CPU vendor")
    set(BUILD_CPU_BRAND    "${OUTPUT_CPU_BRAND}"    CACHE INTERNAL "Build CPU brand")
    set(BUILD_CPU_FAMILY   "${OUTPUT_CPU_FAMILY}"   CACHE INTERNAL "Build CPU family")
    set(BUILD_CPU_MODEL    "${OUTPUT_CPU_MODEL}"    CACHE INTERNAL "Build CPU model")
    set(BUILD_CPU_STEPPING "${OUTPUT_CPU_STEPPING}" CACHE INTERNAL "Build CPU stepping")
    set(BUILD_CPU_FEATURES "${OUTPUT_CPU_FEATURES}" CACHE INTERNAL "Build CPU features")

    ENDIF(NOT DEFINED BUILD_USER)
endmacro(gmx_set_build_information)
