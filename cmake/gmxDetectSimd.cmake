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
            if(CPU_DETECTION_FEATURES MATCHES " avx512er ")
                set(OUTPUT_SIMD "AVX_512_KNL")
            elseif(CPU_DETECTION_FEATURES MATCHES " avx512f ")
                gmx_detect_avx_512_fma_units(NUMBER_OF_AVX_512_FMA_UNITS)
                if(NUMBER_OF_AVX_512_FMA_UNITS EQUAL 2)
                    set(OUTPUT_SIMD "AVX_512")
                elseif(NUMBER_OF_AVX_512_FMA_UNITS EQUAL 1)
                    if (NOT SUGGEST_SIMD_QUIETLY)
                        message(STATUS "This host supports AVX-512, but only has 1 AVX-512 FMA unit, so AVX2 will be faster.")
                    endif()
                    set(OUTPUT_SIMD "AVX2_256")
                else()
                    if (NOT SUGGEST_SIMD_QUIETLY)
                        message(STATUS "Could not run code to detect number of AVX-512 FMA units - assuming 2.")
                    endif()
                    set(OUTPUT_SIMD "AVX_512")
                endif()
            elseif(CPU_DETECTION_FEATURES MATCHES " avx2 ")
                if(CPU_DETECTION_FEATURES MATCHES " amd ")
                    set(OUTPUT_SIMD "AVX2_128")
                else()
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
            elseif(CPU_DETECTION_FEATURES MATCHES " vmx ")
                set(OUTPUT_SIMD "IBM_VMX")
            elseif(CPU_DETECTION_FEATURES MATCHES " qpx ")
                set(OUTPUT_SIMD "IBM_QPX")
            elseif(CPU_DETECTION_FEATURES MATCHES " neon_asimd ")
                set(OUTPUT_SIMD "ARM_NEON_ASIMD")
            elseif(CPU_DETECTION_FEATURES MATCHES " neon " AND NOT GMX_DOUBLE)
                set(OUTPUT_SIMD "ARM_NEON")
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
        if(GMX_TARGET_BGQ)
            # BG/Q requires cross-compilation, so needs this
            # logic. While the qpx feature flag in cpuinfo works, it
            # can't be returned by cpuinfo running on the build host.
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

        string(TOUPPER "${${_suggested_simd}}" ${_suggested_simd})
        set(${_suggested_simd} ${${_suggested_simd}} PARENT_SCOPE)
    endif()
endfunction()
