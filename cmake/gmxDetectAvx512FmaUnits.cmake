#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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
include(gmxSimdFlags)

# gmx_detect_avx_512_fma_units()
#
# Try to detect whether the host has one or two AVX-512 FMA units
# by executing a small program. This will only work on hosts that
# support AVX-512. If successful it sets RESULT to 1 or 2 for the
# number of AVX-512 FMA units, and otherwise -1.
#
function(gmx_detect_avx_512_fma_units RESULT)
    if(CMAKE_CROSSCOMPILING)
        set(${RESULT} -1 CACHE INTERNAL "Result of test for number of AVX-512 FMA units")
    else()
        set(AVX_512_FMA_UNIT_DETECTION_BINARY "${PROJECT_BINARY_DIR}/CMakeFiles/GmxDetectAvx512FmaUnits${CMAKE_EXECUTABLE_SUFFIX}")
        if(NOT AVX_512_FMA_UNIT_DETECTION_COMPILED)

            # Find flags required for AVX-512
            gmx_find_simd_avx_512_flags(SIMD_AVX_512_C_SUPPORTED SIMD_AVX_512_CXX_SUPPORTED
                                        SIMD_AVX_512_C_FLAGS SIMD_AVX_512_CXX_FLAGS)
            # Find flag for GCC inline assembly
            gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

            if(SIMD_AVX_512_CXX_SUPPORTED AND GMX_X86_GCC_INLINE_ASM)
                # Compile the detection program

                set(_compile_definitions "-I${PROJECT_SOURCE_DIR}/src -DGMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE -DSIMD_AVX_512_CXX_SUPPORTED=1 -DGMX_X86_GCC_INLINE_ASM=1 ${SIMD_AVX_512_CXX_FLAGS} ${GMX_STDLIB_CXX_FLAGS}")
                try_compile(AVX_512_FMA_UNIT_DETECTION_COMPILED
                    "${PROJECT_BINARY_DIR}"
                    "${PROJECT_SOURCE_DIR}/src/gromacs/hardware/identifyavx512fmaunits.cpp"
                    COMPILE_DEFINITIONS "${_compile_definitions}"
                    LINK_LIBRARIES "${GMX_STDLIB_LIBRARIES}"
                    OUTPUT_VARIABLE AVX_512_FMA_UNIT_DETECTION_COMPILED_OUTPUT
                    COPY_FILE ${AVX_512_FMA_UNIT_DETECTION_BINARY})
                if(NOT AVX_512_FMA_UNIT_DETECTION_COMPILED AND NOT RUN_AVX_512_FMA_UNIT_DETECTION_COMPILATION_QUIETLY)
                  message(STATUS "Could not identify number of AVX-512 units - detection program did not compile")
                endif()
                set(RUN_AVX_512_FMA_UNIT_DETECTION_COMPILATION_QUIETLY TRUE CACHE INTERNAL "Keep quiet on any future compilation attempts")
            else()
                message(STATUS "Could not identify number of AVX-512 units - detection program missing compilation prerequisites")
            endif()

            if(AVX_512_FMA_UNIT_DETECTION_COMPILED)
                # Run the program
                if(NOT DEFINED ${RESULT})
                    execute_process(COMMAND ${AVX_512_FMA_UNIT_DETECTION_BINARY}
                        RESULT_VARIABLE RESULT_VAR
                        OUTPUT_VARIABLE OUTPUT_VAR_TEMP
                        ERROR_QUIET)
                    if (RESULT_VAR EQUAL 0)
                        string(STRIP "${OUTPUT_VAR_TEMP}" OUTPUT_VAR)
                        set(${RESULT} ${OUTPUT_VAR_TEMP} CACHE INTERNAL "Result of test for number of AVX-512 FMA units")
                    else()
                        message(STATUS "Could not identify number of AVX-512 units - detection program did run successfully")
                        set(${RESULT} -1 CACHE INTERNAL "Result of test for number of AVX-512 FMA units")
                    endif()
                endif()
            endif()
        endif()
      endif()
endfunction()
