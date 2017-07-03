#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

include(gmxTestInlineASM)

# Ensure things like GMX_TARGET_X86 are available
include(gmxDetectTargetArchitecture)
gmx_detect_target_architecture()

# gmx_run_cpu_detection()
#
# Try to detect information about the CPU of the build host by
# building and running the same detection code used by mdrun. This
# works on all architectures where we are not cross-compiling;
# depending on the architecture the detection will either use special
# assembly instructions (like cpuid), preprocessor defines, or probing
# /proc/cpuinfo on Linux.
#
# The TYPE argument is passed as a command-line argument to the
# detection program, and the terminal output is captured and stored in
# the cache variable BUILD_CPU_${TYPE} as the result.
#
# Status messages are written, unless QUIETLY is true.
#
# The function caches information about whether the detection program
# has already been build or run with this TYPE, so this function
# should be called freely, even if the call might be repeated within
# or across invocations of cmake.
#
function(gmx_run_cpu_detection TYPE QUIETLY)
    set(OUTPUT_VAR "")
    # We need to execute the binary, so this only works if not
    # cross-compiling. However, note that we are NOT limited to x86.
    if(CMAKE_CROSSCOMPILING)
        set(EXPLANATION "because we are cross compiling")
    else()
        set(BUILD_CPU_BINARY "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/GmxDetectCpu${CMAKE_EXECUTABLE_SUFFIX}")
        if(NOT BUILD_CPU_COMPILED)
            # Compile the detection program
            set(GMX_TARGET_X86_VALUE 0)
            if(GMX_TARGET_X86)
                set(GMX_TARGET_X86_VALUE 1)
            endif()

            # for x86 we need inline asm to use cpuid
            gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)
            if(GMX_X86_GCC_INLINE_ASM)
                set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM=1")
            else()
                set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM=0")
            endif()

            set(_compile_definitions "${GCC_INLINE_ASM_DEFINE} -I${CMAKE_SOURCE_DIR}/src -DGMX_CPUINFO_STANDALONE ${GMX_STDLIB_CXX_FLAGS} -DGMX_TARGET_X86=${GMX_TARGET_X86_VALUE}")
            set(LINK_LIBRARIES "${GMX_STDLIB_LIBRARIES}")
            try_compile(BUILD_CPU_COMPILED
                "${CMAKE_CURRENT_BINARY_DIR}"
                "${CMAKE_CURRENT_SOURCE_DIR}/src/gromacs/hardware/cpuinfo.cpp"
                COMPILE_DEFINITIONS "${_compile_definitions}"
                CMAKE_FLAGS "-DLINK_LIBRARIES=${LINK_LIBRARIES}"
                OUTPUT_VARIABLE BUILD_CPU_COMPILED_OUTPUT
                COPY_FILE ${BUILD_CPU_BINARY})
        endif()

        if(NOT BUILD_CPU_COMPILED)
            set(EXPLANATION "detection program did not compile")
        else()
            # Run the detection program with -type as the argument.
            string(TOUPPER ${TYPE} UPPERTYPE)
            string(TOLOWER ${TYPE} LOWERTYPE)

            if(NOT DEFINED BUILD_CPU_RUN_${UPPERTYPE})
                execute_process(COMMAND ${BUILD_CPU_BINARY} "-${LOWERTYPE}"
                    RESULT_VARIABLE RESULT_VAR
                    OUTPUT_VARIABLE OUTPUT_VAR_TEMP
                    ERROR_QUIET)
                if (RESULT_VAR EQUAL 0)
                    string(STRIP "${OUTPUT_VAR_TEMP}" OUTPUT_VAR)
                else()
                    set(EXPLANATION "detection program did not run successfully")
                endif()
            endif()
        endif()
    endif()
    set(BUILD_CPU_${UPPERTYPE} "${OUTPUT_VAR}" CACHE INTERNAL "Result of running cpu detection code with argument -${LOWERTYPE}")
    if(NOT ${QUIETLY})
        if (OUTPUT_VAR)
            message(STATUS "Detected build CPU ${TYPE} - ${OUTPUT_VAR}")
        else()
            message(STATUS "Did not detect build CPU ${TYPE} - ${EXPLANATION}")
        endif()
    endif()
endfunction()
