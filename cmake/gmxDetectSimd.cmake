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
    if(${_suggested_simd})
        # There's already been a suggestion made, which can't change
        return()
    endif()

    # for x86 we need inline asm to use cpuid
    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM=1")
    else()
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM=0")
    endif()

    message(STATUS "Detecting best SIMD instructions for this CPU")

    # Get CPU SIMD properties information
    set(_compile_definitions "${GCC_INLINE_ASM_DEFINE} -I${CMAKE_SOURCE_DIR}/src -DGMX_CPUINFO_STANDALONE ${GMX_STDLIB_CXX_FLAGS}")
    if(GMX_TARGET_X86)
        set(_compile_definitions "${_compile_definitions} -DGMX_TARGET_X86")
    endif()

    # Prepare a default suggestion
    set(OUTPUT_SIMD "None")

    # We need to execute the binary, so this only works if not cross-compiling.
    # However, note that we are NOT limited to x86.
    if(NOT CMAKE_CROSSCOMPILING)
        # TODO Extract this try_compile to a helper function, because
        # it duplicates code in gmxSetBuildInformation.cmake
        set(GMX_DETECTSIMD_BINARY "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/GmxDetectSimd${CMAKE_EXECUTABLE_SUFFIX}")
        set(LINK_LIBRARIES "${GMX_STDLIB_LIBRARIES}")
        try_compile(GMX_DETECTSIMD_COMPILED
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/gromacs/hardware/cpuinfo.cpp"
            COMPILE_DEFINITIONS "${_compile_definitions}"
            CMAKE_FLAGS "-DLINK_LIBRARIES=${LINK_LIBRARIES}"
            OUTPUT_VARIABLE GMX_DETECTSIMD_COMPILED_OUTPUT
            COPY_FILE ${GMX_DETECTSIMD_BINARY})
        unset(_compile_definitions)

        if(GMX_DETECTSIMD_COMPILED)
            if(NOT DEFINED GMX_DETECTSIMD_RUN)
                execute_process(COMMAND ${GMX_DETECTSIMD_BINARY} "-features"
                    RESULT_VARIABLE GMX_DETECTSIMD_RUN
                    OUTPUT_VARIABLE OUTPUT_TMP
                    ERROR_QUIET)
                set(GMX_DETECTSIMD_RUN "${GMX_DETECTSIMD_RUN}" CACHE INTERNAL "Result of running cpuinfo code to detect SIMD support")
                if(GMX_DETECTSIMD_RUN EQUAL 0)
                    # Make a concrete suggestion of SIMD level
                    if(GMX_TARGET_X86)
                        if(OUTPUT_TMP MATCHES " avx512er ")
                            set(OUTPUT_SIMD "AVX_512_KNL")
                        elseif(OUTPUT_TMP MATCHES " avx512f ")
                            set(OUTPUT_SIMD "AVX_512")
                        elseif(OUTPUT_TMP MATCHES " avx2 ")
                            set(OUTPUT_SIMD "AVX2_256")
                        elseif(OUTPUT_TMP MATCHES " avx ")
                            if(OUTPUT_TMP MATCHES " fma4 ")
                                # AMD that works better with avx-128-fma
                                set(OUTPUT_SIMD "AVX_128_FMA")
                            else()
                                # Intel
                                set(OUTPUT_SIMD "AVX_256")
                            endif()
                        elseif(OUTPUT_TMP MATCHES " sse4.1 ")
                            set(OUTPUT_SIMD "SSE4.1")
                        elseif(OUTPUT_TMP MATCHES " sse2 ")
                            set(OUTPUT_SIMD "SSE2")
                        endif()
                    else()
                        if(OUTPUT_TMP MATCHES " vsx ")
                            set(OUTPUT_SIMD "IBM_VSX")
                        elseif(OUTPUT_TMP MATCHES " vmx ")
                            set(OUTPUT_SIMD "IBM_VMX")
                        elseif(OUTPUT_TMP MATCHES " qpx ")
                            set(OUTPUT_SIMD "IBM_QPX")
                        elseif(OUTPUT_TMP MATCHES " neon_asimd ")
                            set(OUTPUT_SIMD "ARM_NEON_ASIMD")
                        elseif(OUTPUT_TMP MATCHES " neon ")
                            set(OUTPUT_SIMD "ARM_NEON")
                        endif()
                    endif()
                    message(STATUS "Detected best SIMD instructions for this CPU - ${OUTPUT_SIMD}")
                else()
                    message(WARNING "Cannot run cpuinfo code, which means no SIMD suggestion can be made.")
                    message(STATUS "Run output: ${OUTPUT_TMP}")
                endif()
            endif()
        else()
            message(WARNING "Cannot compile cpuinfo code, which means no SIMD instructions.")
            message(STATUS "Compile output: ${GMX_DETECTSIMD_COMPILED_OUTPUT}")
        endif()
    else()
        message(WARNING "Cannot detect SIMD architecture for this cross-compile; you should check it manually.")
    endif()

    set(${_suggested_simd} "${OUTPUT_SIMD}" CACHE INTERNAL "Suggested SIMD")
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
