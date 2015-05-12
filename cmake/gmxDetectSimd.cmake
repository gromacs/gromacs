#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
# gmx_detect_simd(GMX_SUGGESTED_SIMD)
#
# Try to detect CPU information and suggest SIMD instruction set
# that fits the current CPU. This should work on all architectures
# where we are not cross-compiling; depending on the architecture the
# detection will either use special assembly instructions (like cpuid),
# preprocessor defines, or probing /proc/cpuinfo on Linux.
# 
# This assumes gmx_detect_target_architecture() has already been run,
# so that things like GMX_TARGET_X86 are already available.
# (otherwise we cannot use inline ASM on x86).
#
# Sets ${GMX_SUGGESTED_SIMD} in the parent scope if
# GMX_SIMD is not set (e.g. by the user, or a previous run
# of CMake).
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)

function(gmx_suggest_simd _suggested_simd)

    # for x86 we need inline asm to use cpuid
    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM")
    else()
        set(GCC_INLINE_ASM_DEFINE "")
    endif()

    message(STATUS "Detecting best SIMD instructions for this CPU")

    # Get CPU SIMD properties information
    set(_compile_definitions "${GCC_INLINE_ASM_DEFINE} -I${CMAKE_SOURCE_DIR}/src -DGMX_CPUID_STANDALONE")
    if(GMX_TARGET_X86)
        set(_compile_definitions "${_compile_definitions} -DGMX_TARGET_X86")
    endif()

    # We need to execute the binary, so this only works if not cross-compiling.
    # However, note that we are NOT limited to x86.
    if(NOT CMAKE_CROSSCOMPILING)
        try_run(GMX_CPUID_RUN_SIMD GMX_CPUID_COMPILED
                ${CMAKE_BINARY_DIR}
                ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
                COMPILE_DEFINITIONS ${_compile_definitions}
                RUN_OUTPUT_VARIABLE OUTPUT_TMP
                COMPILE_OUTPUT_VARIABLE GMX_CPUID_COMPILE_OUTPUT
                ARGS "-simd")

        if(NOT GMX_CPUID_COMPILED)
            message(WARNING "Cannot compile CPUID code, which means no SIMD instructions.")
            message(STATUS "Compile output: ${GMX_CPUID_COMPILE_OUTPUT}")
            set(OUTPUT_TMP "None")
        elseif(NOT GMX_CPUID_RUN_SIMD EQUAL 0)
            message(WARNING "Cannot run CPUID code, which means no SIMD instructions.")
            message(STATUS "Run output: ${OUTPUT_TMP}")
            set(OUTPUT_TMP "None")
        endif(NOT GMX_CPUID_COMPILED)

        string(STRIP "${OUTPUT_TMP}" OUTPUT_SIMD)

        set(${_suggested_simd} "${OUTPUT_SIMD}" PARENT_SCOPE)
        message(STATUS "Detected best SIMD instructions for this CPU - ${OUTPUT_SIMD}")
    else()
        set(${_suggested_simd} "None" PARENT_SCOPE)
        message(WARNING "Cannot detect SIMD architecture for this cross-compile; you should check it manually.")
    endif()
endfunction()

function(gmx_detect_simd _suggested_simd)
    if(NOT DEFINED GMX_SIMD)
        if(GMX_TARGET_BGQ)
            set(${_suggested_simd} "IBM_QPX")
        elseif(GMX_TARGET_FUJITSU_SPARC64)
            # HPC-ACE is always present. In the future we
            # should add detection for HPC-ACE2 here.
            set(${_suggested_simd} "Sparc64_HPC_ACE")
        elseif(GMX_TARGET_MIC)
            set(${_suggested_simd} "MIC")
        else()
            gmx_suggest_simd(${_suggested_simd})
        endif()

        set(${_suggested_simd} ${${_suggested_simd}} PARENT_SCOPE)
    endif()
endfunction()
