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
# - Check the username performing the build, as well as date and time
#
# gmx_detect_acceleration(GMX_SUGGESTED_CPU_ACCELERATION)
#
# Try to detect CPU information and suggest an acceleration option
# (such as SSE/AVX) that fits the current CPU. These functions assume
# that gmx_detect_target_architecture() has already been run, so that
# things like GMX_IS_X86 are already available.
#
# Sets ${GMX_SUGGESTED_CPU_ACCELERATION} in the parent scope if
# GMX_CPU_ACCELERATION is not set (e.g. by the user, or a previous run
# of CMake).
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)

function(gmx_suggest_x86_acceleration _suggested_acceleration)

    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM")
    else(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "")
    endif(GMX_X86_GCC_INLINE_ASM)

    message(STATUS "Detecting best acceleration for this CPU")

    # Get CPU acceleration information
    try_run(GMX_CPUID_RUN_ACC GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE -DGMX_IS_X86"
            RUN_OUTPUT_VARIABLE OUTPUT_TMP
            COMPILE_OUTPUT_VARIABLE GMX_CPUID_COMPILE_OUTPUT
            ARGS "-acceleration")

    if(NOT GMX_CPUID_COMPILED)
        message(WARNING "Cannot compile CPUID code, which means no CPU-specific acceleration.")
        message(STATUS "Compile output: ${GMX_CPUID_COMPILE_OUTPUT}")
        set(OUTPUT_TMP "None")
    elseif(NOT GMX_CPUID_RUN_ACC EQUAL 0)
        message(WARNING "Cannot run CPUID code, which means no CPU-specific optimization.")
        message(STATUS "Run output: ${OUTPUT_TMP}")
        set(OUTPUT_TMP "None")
    endif(NOT GMX_CPUID_COMPILED)

    string(STRIP "@OUTPUT_TMP@" OUTPUT_ACC)

    set(${_suggested_acceleration} "@OUTPUT_ACC@" PARENT_SCOPE)
    message(STATUS "Detected best acceleration for this CPU - @OUTPUT_ACC@")
endfunction()

function(gmx_detect_acceleration _suggested_acceleration)
    if(NOT DEFINED GMX_CPU_ACCELERATION)
        if(GMX_IS_BGQ)
            set(${_suggested_acceleration} "IBM_QPX")
        elseif(GMX_IS_X86)
            gmx_suggest_x86_acceleration(${_suggested_acceleration})
        else()
            set(${_suggested_acceleration} "None")
        endif()

        set(${_suggested_acceleration} ${${_suggested_acceleration}} PARENT_SCOPE)
    endif()
endfunction()
