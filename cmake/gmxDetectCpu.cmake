#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017,2019, by the GROMACS development team, led by
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
# the cache variable CPU_DETECTION_${TYPE} as the result. If the detection
# program fails to compile, or fails to run, no value is stored.
#
# The function caches information about whether the detection program
# has already been built or run with this TYPE, so this function
# should be called freely, even if the call might be repeated within
# or across invocations of cmake.
#
function(gmx_run_cpu_detection TYPE)
    string(TOUPPER ${TYPE} UPPERTYPE)
    string(TOLOWER ${TYPE} LOWERTYPE)

    set(OUTPUT_VAR "")
    # We need to execute the binary, so this only works if not
    # cross-compiling. However, note that we are NOT limited to x86.
    if(CMAKE_CROSSCOMPILING)
        # TODO Need we explain that we're not detecting because we are cross compiling?
    else()
        set(CPU_DETECTION_BINARY "${PROJECT_BINARY_DIR}/CMakeFiles/GmxDetectCpu${CMAKE_EXECUTABLE_SUFFIX}")
        if(NOT CPU_DETECTION_COMPILED)
            # Compile the detection program
            set(GMX_TARGET_X86_VALUE 0)
            if(GMX_TARGET_X86)
                set(GMX_TARGET_X86_VALUE 1)
            endif()

            # for x86 we need inline assembly to use cpuid
            gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)
            if(GMX_X86_GCC_INLINE_ASM)
                set(GCC_INLINE_ASM_DEFINE -DGMX_X86_GCC_INLINE_ASM=1)
            else()
                set(GCC_INLINE_ASM_DEFINE -DGMX_X86_GCC_INLINE_ASM=0)
            endif()

            set(_compile_definitions ${GCC_INLINE_ASM_DEFINE};-I${PROJECT_SOURCE_DIR}/src;-DGMX_CPUINFO_STANDALONE=1;-DGMX_TARGET_X86=${GMX_TARGET_X86_VALUE})
            try_compile(CPU_DETECTION_COMPILED
                "${PROJECT_BINARY_DIR}"
                "${PROJECT_SOURCE_DIR}/src/gromacs/hardware/cpuinfo.cpp"
                COMPILE_DEFINITIONS "${_compile_definitions}"
                CMAKE_FLAGS "-DLINK_LIBRARIES=${LINK_LIBRARIES}"
                OUTPUT_VARIABLE CPU_DETECTION_COMPILED_OUTPUT
                COPY_FILE ${CPU_DETECTION_BINARY})
            if(NOT CPU_DETECTION_COMPILED AND NOT RUN_CPU_DETECTION_COMPILATION_QUIETLY)
                if(GMX_TARGET_X86)
                    message(WARNING "CPU detection program did not compile on x86 host - this should never happen. It is VERY bad for performance, since you will lose all SIMD support. Please file a bug report.")
                else()
                    message(WARNING "Did not detect build CPU ${LOWERTYPE} - detection program did not compile. Please file a bug report if this is a common platform.")
                endif()
            endif()
            set(RUN_CPU_DETECTION_COMPILATION_QUIETLY TRUE CACHE INTERNAL "Keep quiet on any future compilation attempts")
        endif()

        if(CPU_DETECTION_COMPILED)
            # Run the detection program with -type as the argument.

            if(NOT DEFINED CPU_DETECTION_${UPPERTYPE})
                execute_process(COMMAND ${CPU_DETECTION_BINARY} "-${LOWERTYPE}"
                    RESULT_VARIABLE RESULT_VAR
                    OUTPUT_VARIABLE OUTPUT_VAR_TEMP
                    ERROR_QUIET)
                if (RESULT_VAR EQUAL 0)
                    string(STRIP "${OUTPUT_VAR_TEMP}" OUTPUT_VAR)
                    message(STATUS "Detected build CPU ${LOWERTYPE} - ${OUTPUT_VAR}")
                    set(CPU_DETECTION_${UPPERTYPE} "${OUTPUT_VAR}" CACHE INTERNAL "Result of running cpu detection code with argument -${LOWERTYPE}")
                else()
                    message(STATUS "Did not detect build CPU ${LOWERTYPE} - detection program did not run successfully")
                endif()
            endif()
        endif()
    endif()
endfunction()
