#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# - Check the username performing the build, as well as date and time
#
# gmx_detect_simd(_suggested_simd)
#
# Try to detect CPU features and suggest a SIMD instruction set
# that fits the current CPU. This should work on all architectures
# where we are not cross-compiling; depending on the architecture the
# detection will either use special assembly instructions (like cpuid),
# preprocessor defines, or probing /proc/cpuinfo on Linux.
# 
# Sets ${_suggested_simd} in the parent scope if GMX_SIMD is not set
# (e.g. by the user, or a previous run of CMake).
# The string is converted to uppercase for compatibility with
# gmx_option_multichoice() user input parsing.
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)
# Ensure things like GMX_TARGET_X86 are available
include(gmxDetectTargetArchitecture)
gmx_detect_target_architecture()

include(gmxDetectCpu)
include(gmxDetectAvx512FmaUnits)

function(gmx_suggest_simd _suggested_simd)
    if (NOT SUGGEST_SIMD_QUIETLY)
        message(STATUS "Detecting best SIMD instructions for this CPU")
    endif()

    # Prepare a default suggestion
    set(OUTPUT_SIMD "None")

    # Detect CPU features and place the string in CPU_DETECTION_FEATURES
    # Note that we are NOT limited to x86.
    gmx_run_cpu_detection(features)

    if (DEFINED CPU_DETECTION_FEATURES)
        # Make a concrete suggestion of SIMD level if a feature flag
        # matches. Make sure that the match strings below work even if
        # the feature is first or last.
        set(CPU_DETECTION_FEATURES " ${CPU_DETECTION_FEATURES} ")

        if(GMX_TARGET_X86)
            gmx_run_cpu_detection(brand)
            if(CPU_DETECTION_FEATURES MATCHES " avx512er ")
                set(OUTPUT_SIMD "AVX_512_KNL")
            elseif(CPU_DETECTION_FEATURES MATCHES " avx512f ")
                if(CPU_DETECTION_BRAND MATCHES "Intel")
                    gmx_detect_avx_512_fma_units(NUMBER_OF_AVX_512_FMA_UNITS)
                    if(NUMBER_OF_AVX_512_FMA_UNITS EQUAL 2)
                        set(OUTPUT_SIMD "AVX_512")
                    elseif(NUMBER_OF_AVX_512_FMA_UNITS EQUAL 1)
                        if (NOT SUGGEST_SIMD_QUIETLY)
                            message(STATUS "This is an Intel CPU with only 1 AVX-512 FMA unit, so AVX2 will be faster.")
                        endif()
                        set(OUTPUT_SIMD "AVX2_256")
                    else()
                        if (NOT SUGGEST_SIMD_QUIETLY)
                            message(STATUS "Could not run code to detect number of AVX-512 FMA units - assuming 2.")
                        endif()
                        set(OUTPUT_SIMD "AVX_512")
                    endif()
                else()
                    # Non-Intel vendor with AVX-512 presently means AMD,
                    # and this far AVX-512 is always faster on AMD, even with a single FMA unit.
                    set(OUTPUT_SIMD "AVX_512")
                endif()
            elseif(CPU_DETECTION_FEATURES MATCHES " avx2 ")
                if(CPU_DETECTION_BRAND MATCHES "AMD")
                    gmx_run_cpu_detection(family)
                    gmx_run_cpu_detection(model)
                    set(ZEN1_MODELS 1 17 8 24)
                    if("${CPU_DETECTION_FAMILY}" STREQUAL "23" AND "${CPU_DETECTION_MODEL}" IN_LIST ZEN1_MODELS)
                        # Zen/Zen+, where 128-bit AVX2 will be faster
                        set(OUTPUT_SIMD "AVX2_128")
                    else()
                        # Zen2 or later, where 256-bit AVX2 should be faster
                        set(OUTPUT_SIMD "AVX2_256")
                    endif()
                else()
                    # not AMD
                    set(OUTPUT_SIMD "AVX2_256")
                endif()
            elseif(CPU_DETECTION_FEATURES MATCHES " avx ")
                if(CPU_DETECTION_FEATURES MATCHES " fma4 ")
                    # AMD that works better with avx-128-fma
                    set(OUTPUT_SIMD "AVX_128_FMA")
                else()
                    # Intel
                    set(OUTPUT_SIMD "AVX_256")
                endif()
            elseif(CPU_DETECTION_FEATURES MATCHES " sse4.1 ")
                set(OUTPUT_SIMD "SSE4.1")
            elseif(CPU_DETECTION_FEATURES MATCHES " sse2 ")
                set(OUTPUT_SIMD "SSE2")
            endif()
        else()
            if(CPU_DETECTION_FEATURES MATCHES " vsx ")
                set(OUTPUT_SIMD "IBM_VSX")
            elseif(CPU_DETECTION_FEATURES MATCHES " sve ")
                set(OUTPUT_SIMD "ARM_SVE")
            elseif(CPU_DETECTION_FEATURES MATCHES " neon_asimd ")
                set(OUTPUT_SIMD "ARM_NEON_ASIMD")
            endif()
        endif()
        if (NOT SUGGEST_SIMD_QUIETLY)
            message(STATUS "Detected best SIMD instructions for this CPU - ${OUTPUT_SIMD}")
        endif()
    else()
        if (NOT SUGGEST_SIMD_QUIETLY)
            message(STATUS "Detection for best SIMD instructions failed, using SIMD - ${OUTPUT_SIMD}")
        endif()
    endif()

    set(${_suggested_simd} "${OUTPUT_SIMD}" PARENT_SCOPE)
    set(SUGGEST_SIMD_QUIETLY TRUE CACHE INTERNAL "Be quiet during future construction of SIMD suggestions")
endfunction()

function(gmx_detect_simd _suggested_simd)
    if(GMX_SIMD STREQUAL "AUTO")
        gmx_suggest_simd(${_suggested_simd})

        string(TOUPPER "${${_suggested_simd}}" ${_suggested_simd})
        set(${_suggested_simd} ${${_suggested_simd}} PARENT_SCOPE)
    endif()
endfunction()
