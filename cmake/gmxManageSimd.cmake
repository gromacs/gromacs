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

include(gmxDetectCpu)
include(gmxFindFlagsForSource)

# Macro that manages setting the respective C and C++ toolchain
# variables so that subsequent tests for SIMD support can work.
macro(prepare_x86_toolchain TOOLCHAIN_C_FLAGS_VARIABLE TOOLCHAIN_CXX_FLAGS_VARIABLE)
    # On OS X, we often want to use gcc instead of clang, since gcc
    # supports OpenMP (until clang 3.8, or so, plus whenever Apple
    # support it in their version). However, by default gcc uses the
    # external system assembler, which does not support AVX, so we
    # need to tell the linker to use the clang compilers assembler
    # instead - and this has to happen before we detect AVX flags.
    if(APPLE AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
        gmx_test_cflag(GNU_C_USE_CLANG_AS "-Wa,-q" ${TOOLCHAIN_C_FLAGS_VARIABLE})
    endif()
    if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        gmx_test_cxxflag(GNU_CXX_USE_CLANG_AS "-Wa,-q" ${TOOLCHAIN_CXX_FLAGS_VARIABLE})
    endif()
endmacro()

# Macro that manages setting the respective C and C++ toolchain
# variables so that subsequent tests for SIMD support can work.
macro(prepare_power_vsx_toolchain TOOLCHAIN_C_FLAGS_VARIABLE TOOLCHAIN_CXX_FLAGS_VARIABLE)
    if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_C_COMPILER_ID} MATCHES "GNU")
        # VSX uses the same function API as Altivec/VMX, so make sure we tune for the current CPU and not VMX.
        # By putting these flags here rather than in the general compiler flags file we can safely assume
        # that we are at least on Power7 since that is when VSX appeared.
        gmx_run_cpu_detection(brand)
        if(CPU_DETECTION_BRAND MATCHES "POWER7")
            gmx_test_cflag(GNU_C_VSX_POWER7   "-mcpu=power7 -mtune=power7" ${TOOLCHAIN_C_FLAGS_VARIABLE})
            gmx_test_cflag(GNU_CXX_VSX_POWER7 "-mcpu=power7 -mtune=power7" ${TOOLCHAIN_CXX_FLAGS_VARIABLE})
        else()
            # Enable power8 vector extensions on all platforms except old Power7.
            gmx_test_cflag(GNU_C_VSX_POWER8   "-mcpu=power8 -mpower8-vector -mpower8-fusion -mdirect-move" ${TOOLCHAIN_C_FLAGS_VARIABLE})
            gmx_test_cflag(GNU_CXX_VSX_POWER8 "-mcpu=power8 -mpower8-vector -mpower8-fusion -mdirect-move" ${TOOLCHAIN_CXX_FLAGS_VARIABLE})
        endif()
        # Altivec was originally single-only, and it took a while for compilers
        # to support the double-precision features in VSX.
        if(GMX_DOUBLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
            message(FATAL_ERROR "Using VSX SIMD in double precision with GCC requires GCC-4.9 or later.")
        endif()
    endif()
    if(${CMAKE_CXX_COMPILER_ID} MATCHES "XL" OR ${CMAKE_C_COMPILER_ID} MATCHES "XL")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13.1.5" OR CMAKE_C_COMPILER_VERSION VERSION_LESS "13.1.5")
            message(FATAL_ERROR "Using VSX SIMD requires XL compiler version 13.1.5 or later.")
        endif()
    endif()
endmacro()

# Issue a fatal error with an appropriate message, when the toolchain
# was not able to compile code for SIMD support.
#
# Inputs:
#  SIMD_STRING              A string describing the kind of SIMD support that didn't work.
#  ALTERNATIVE_SUGGESTION   A string describing anything the user could try other than getting a new compiler.
#  SUGGEST_BINUTILS_UPDATE  True when there's information that the compiler was OK, but something else was not.
function(gmx_give_fatal_error_when_simd_support_not_found SIMD_STRING ALTERNATIVE_SUGGESTION SUGGEST_BINUTILS_UPDATE)
    if(SUGGEST_BINUTILS_UPDATE)
        set(_msg "Found a compiler flag for ${SIMD_STRING} support, but some other problem exists. Update your assembler and/or linker, e.g. in the binutils package of your distribution.")
    else()
        set(_msg "Cannot find ${SIMD_STRING} compiler flag. Use a newer compiler, or ${ALTERNATIVE_SUGGESTION}.")
    endif()
    message(FATAL_ERROR ${_msg})
endfunction()

macro(gmx_manage_simd)

set(GMX_SIMD_ACCURACY_BITS_SINGLE 22 CACHE STRING "Target mantissa bits for SIMD single math")
#
# Note that we typically restrict double precision target accuracy to be twice that
# of single. This means we only need one more N-R iteration for 1/sqrt(x) and 1(x),
# and the first iteration can sometimes be done as a pair in single precision. This should
# be plenty enough for Molecular Dynamics applications. Many of our double precision math
# functions still achieve very close to full double precision, but we do not guarantee that
# they will be able to achieve higher accuracy if you set this beyond 44 bits. GROMACS will
# work - but some unit tests might fail.
#
set(GMX_SIMD_ACCURACY_BITS_DOUBLE 44 CACHE STRING "Target mantissa bits for SIMD double math")
mark_as_advanced(GMX_SIMD_ACCURACY_BITS_SINGLE)
mark_as_advanced(GMX_SIMD_ACCURACY_BITS_DOUBLE)

if(${GMX_SIMD_ACCURACY_BITS_SINGLE} GREATER 22)
    message(STATUS "Note: Full mantissa accuracy (including least significant bit) requested for SIMD single math. Presently we cannot get the least significant bit correct since that would require different algorithms - reducing to 22 bits.")
    set(GMX_SIMD_ACCURACY_BITS_SINGLE 22 CACHE STRING "Target mantissa bits for SIMD single math" FORCE)
endif()

if(${GMX_SIMD_ACCURACY_BITS_DOUBLE} GREATER 51)
    message(STATUS "Note: Full mantissa accuracy (including least significant bit) requested for SIMD double math. Presently we cannot get the least significant bit correct since that would require different algorithms - reducing to 51 bits.")
    set(GMX_SIMD_ACCURACY_BITS_DOUBLE 51 CACHE STRING "Target mantissa bits for SIMD double math" FORCE)
endif()

#
# Section to set (and test) compiler flags for SIMD.
#
# If the user chose the (default) automatic behaviour, then detection
# is run to suggest a SIMD choice suitable for the build
# host. Otherwise, the users's choice is always honoured. The compiler
# flags will be set based on that choice.
#

set(GMX_SIMD_ACTIVE ${GMX_SIMD})
if(GMX_SIMD STREQUAL "AUTO")
    include(gmxDetectSimd)
    gmx_detect_simd(GMX_SUGGESTED_SIMD)
    set(GMX_SIMD_ACTIVE ${GMX_SUGGESTED_SIMD})
endif()

if(GMX_SIMD_ACTIVE STREQUAL "NONE")
    # nothing to do configuration-wise
    set(SIMD_STATUS_MESSAGE "SIMD instructions disabled")
elseif(GMX_SIMD_ACTIVE STREQUAL "SSE2")

    gmx_find_flags(
        "#include<xmmintrin.h>
         int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_rsqrt_ps(x);return _mm_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-msse2" "/arch:SSE2" "-hgnu")

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("SSE2" "disable SIMD support (slow)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling SSE2 SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "SSE4.1")

    # Note: MSVC enables SSE4.1 with the SSE2 flag, so we include that in testing.
    gmx_find_flags(
        "#include<smmintrin.h>
        int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_dp_ps(x,x,0x77);return _mm_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_SSE_4_1_C_FLAGS SIMD_SSE_4_1_CXX_FLAGS
        "-msse4.1" "/arch:SSE4.1" "/arch:SSE2" "-hgnu")

    if(NOT SIMD_SSE_4_1_C_FLAGS OR NOT SIMD_SSE_4_1_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("SSE4.1" "choose SSE2 SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_SSE4_1 1)
    set(SIMD_STATUS_MESSAGE "Enabling SSE4.1 SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_128_FMA")

    prepare_x86_toolchain(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    # We don't have the full compiler version string yet (BUILD_C_COMPILER),
    # so we can't distinguish vanilla from Apple clang versions, but catering for a few rare AMD
    # hackintoshes is not worth the effort.
    if (APPLE AND (CMAKE_C_COMPILER_ID STREQUAL "Clang" OR
                CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
        message(WARNING "Due to a known compiler bug, Clang up to version 3.2 (and Apple Clang up to version 4.1) produces incorrect code with AVX_128_FMA SIMD. As we cannot work around this bug on OS X, you will have to select a different compiler or SIMD instruction set.")
    endif()

    # clang <=3.2 contains a bug that causes incorrect code to be generated for the
    # vfmaddps instruction and therefore the bug is triggered with AVX_128_FMA.
    # (see: http://llvm.org/bugs/show_bug.cgi?id=15040).
    # We can work around this by not using the integrated assembler (except on OS X
    # which has an outdated assembler that does not support AVX instructions).
    if (CMAKE_C_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_VERSION VERSION_LESS "3.3")
        # we assume that we have an external assembler that supports AVX
        message(STATUS "Clang ${CMAKE_C_COMPILER_VERSION} detected, enabling FMA bug workaround")
        set(TOOLCHAIN_C_FLAGS "${TOOLCHAIN_C_FLAGS} -no-integrated-as")
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
        # we assume that we have an external assembler that supports AVX
        message(STATUS "Clang ${CMAKE_CXX_COMPILER_VERSION} detected, enabling FMA bug workaround")
        set(TOOLCHAIN_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS} -no-integrated-as")
    endif()

    # AVX128/FMA on AMD is a bit complicated. We need to do detection in three stages:
    # 1) Find the flags required for generic AVX support
    # 2) Find the flags necessary to enable fused-multiply add support
    # 3) Optional: Find a flag to enable the AMD XOP instructions

    ### STAGE 1: Find the generic AVX flag
    gmx_find_flags(
        "#include<immintrin.h>
        int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_permute_ps(x,1);return 0;}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_GENERIC_AVX_C_FLAGS SIMD_GENERIC_AVX_CXX_FLAGS
        "-mavx" "/arch:AVX" "-hgnu")

    ### STAGE 2: Find the fused-multiply add flag.
    # GCC requires x86intrin.h for FMA support. MSVC 2010 requires intrin.h for FMA support.
    check_include_file(x86intrin.h HAVE_X86INTRIN_H ${SIMD_C_FLAGS})
    check_include_file(intrin.h HAVE_INTRIN_H ${SIMD_C_FLAGS})
    if(HAVE_X86INTRIN_H)
        set(INCLUDE_X86INTRIN_H "#include <x86intrin.h>")
    endif()
    if(HAVE_INTRIN_H)
        set(INCLUDE_INTRIN_H "#include <xintrin.h>")
    endif()

    gmx_find_flags(
        "#include<immintrin.h>
        ${INCLUDE_X86INTRIN_H}
        ${INCLUDE_INTRIN_H}
        int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_macc_ps(x,x,x);return _mm_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-mfma4" "-hgnu")

    # We only need to check the last (FMA) test; that will always fail if the generic AVX test failed
    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("128-bit AVX with FMA support" "choose SSE4.1 SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    ### STAGE 3: Optional: Find the XOP instruction flag (No point in yelling if this does not work)
    gmx_find_flags(
        "#include<immintrin.h>
        ${INCLUDE_X86INTRIN_H}
        ${INCLUDE_INTRIN_H}
        int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_frcz_ps(x);return _mm_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_AVX_128_XOP_C_FLAGS SIMD_AVX_128_XOP_CXX_FLAGS
        "-mxop")

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 128-bit AVX SIMD GROMACS SIMD (with fused-multiply add)")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_256")

    prepare_x86_toolchain(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(
        "#include<immintrin.h>
         int main(){__m256 x=_mm256_set1_ps(0.5);x=_mm256_add_ps(x,x);return _mm256_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-mavx" "/arch:AVX" "-hgnu")

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("AVX" "choose SSE4.1 SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 256-bit AVX SIMD instructions")

elseif(GMX_SIMD_ACTIVE MATCHES "AVX2_")

    prepare_x86_toolchain(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(
        "#include<immintrin.h>
         int main(){__m256i x=_mm256_set1_epi32(5);x=_mm256_add_epi32(x,x);return _mm256_movemask_epi8(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-march=core-avx2" "-mavx2" "/arch:AVX" "-hgnu") # no AVX2-specific flag for MSVC yet

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("AVX2" "choose AVX SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)

    if(GMX_SIMD_ACTIVE STREQUAL "AVX2_128")
        set(SIMD_STATUS_MESSAGE "Enabling 128-bit AVX2 SIMD instructions")
    else()
        set(SIMD_STATUS_MESSAGE "Enabling 256-bit AVX2 SIMD instructions")
    endif()

elseif(GMX_SIMD_ACTIVE STREQUAL "MIC")

    # No flags needed. Not testing.
    set(GMX_SIMD_X86_MIC 1)
    set(SIMD_STATUS_MESSAGE "Enabling MIC (Xeon Phi) SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_512")

    prepare_x86_toolchain(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(
        "#include<immintrin.h>
         int main(){__m512 y,x=_mm512_set1_ps(0.5);y=_mm512_fmadd_ps(x,x,x);return (int)_mm512_cmp_ps_mask(x,y,_CMP_LT_OS);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-xCORE-AVX512 -qopt-zmm-usage=high" "-xCORE-AVX512" "-mavx512f -mfma" "-mavx512f" "/arch:AVX" "-hgnu") # no AVX_512F flags known for MSVC yet. ICC should use ZMM if code anyhow uses ZMM

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("AVX 512F" "choose a lower level of SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 512-bit AVX-512 SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_512_KNL")

    prepare_x86_toolchain(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(
        "#include<immintrin.h>
        int main(){__m512 y,x=_mm512_set1_ps(0.5);y=_mm512_rsqrt28_ps(x);return (int)_mm512_cmp_ps_mask(x,y,_CMP_LT_OS);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-xMIC-AVX512" "-mavx512er -mfma" "-mavx512er" "/arch:AVX" "-hgnu") # no AVX_512ER flags known for MSVC yet

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("AVX 512ER" "choose a lower level of SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 512-bit AVX-512-KNL SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "ARM_NEON")

    if (GMX_DOUBLE)
        message(FATAL_ERROR "ARM_NEON SIMD support is not available for a double precision build because the architecture lacks double-precision support")
    endif()

    gmx_find_flags(
        "#include<arm_neon.h>
         int main(){float32x4_t x=vdupq_n_f32(0.5);x=vmlaq_f32(x,x,x);return vgetq_lane_f32(x,0)>0;}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-mfpu=neon-vfpv4" "-mfpu=neon" "")

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("ARM NEON" "disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 32-bit ARM NEON SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "ARM_NEON_ASIMD")

    gmx_find_flags(
        "#include<arm_neon.h>
         int main(){float64x2_t x=vdupq_n_f64(0.5);x=vfmaq_f64(x,x,x);x=vrndnq_f64(x);return vgetq_lane_f64(x,0)>0;}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "")

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("ARM (AArch64) NEON Advanced SIMD" "particularly gcc version 4.9 or later, or disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling ARM (AArch64) NEON Advanced SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "IBM_QPX")

    try_compile(TEST_QPX ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestQPX.c")

    if (TEST_QPX)
        message(WARNING "IBM QPX SIMD instructions selected. This will work, but SIMD kernels are only available for the Verlet cut-off scheme. The plain C kernels that are used for the group cut-off scheme kernels will be slow, so please consider using the Verlet cut-off scheme.")
        set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
        set(SIMD_STATUS_MESSAGE "Enabling IBM QPX SIMD instructions")

    else()
        gmx_give_fatal_error_when_simd_support_not_found("IBM QPX" "or 'cmake .. -DCMAKE_TOOLCHAIN_FILE=Platform/BlueGeneQ-static-bgclang-CXX' to set up the tool chain" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

elseif(GMX_SIMD_ACTIVE STREQUAL "IBM_VMX")

    gmx_find_flags(
        "#include<altivec.h>
         int main(){vector float x,y=vec_ctf(vec_splat_s32(1),0);x=vec_madd(y,y,y);return vec_all_ge(y,x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-maltivec -mabi=altivec" "-qarch=auto -qaltivec")

    if(NOT SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS OR NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("IBM VMX" "disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling IBM VMX SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "IBM_VSX")

    prepare_power_vsx_toolchain(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(
        "#include<altivec.h>
         int main(){vector double x,y=vec_splats(1.0);x=vec_madd(y,y,y);return vec_all_ge(y,x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_${GMX_SIMD_ACTIVE}_C_FLAGS SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS
        "-mvsx" "-maltivec -mabi=altivec" "-qarch=auto -qaltivec")

    # Usually we check also for the C compiler here, but a C compiler
    # is not required for SIMD support on this platform. cmake through
    # at least version 3.7 cannot pass this check with the C compiler
    # in the latest xlc 13.1.5, but the C++ compiler has different
    # behaviour and is OK. See Redmine #2102.
    if(NOT SIMD_${GMX_SIMD_ACTIVE}_CXX_FLAGS)
        gmx_give_fatal_error_when_simd_support_not_found("IBM VSX" "disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${TOOLCHAIN_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling IBM VSX SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "SPARC64_HPC_ACE")

    # Note that GMX_RELAXED_DOUBLE_PRECISION is enabled by default in the top-level CMakeLists.txt

    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling Sparc64 HPC-ACE SIMD instructions")

elseif(GMX_SIMD_ACTIVE STREQUAL "REFERENCE")

    # NB: This file handles settings for the SIMD module, so in the interest 
    # of proper modularization, please do NOT put any verlet kernel settings in this file.

    if(GMX_SIMD_REF_FLOAT_WIDTH)
        add_definitions(-DGMX_SIMD_REF_FLOAT_WIDTH=${GMX_SIMD_REF_FLOAT_WIDTH})
    endif()
    if(GMX_SIMD_REF_DOUBLE_WIDTH)
      	add_definitions(-DGMX_SIMD_REF_DOUBLE_WIDTH=${GMX_SIMD_REF_DOUBLE_WIDTH})
    endif()

    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling reference (emulated) SIMD instructions.")

else()
    gmx_invalid_option_value(GMX_SIMD_ACTIVE)
endif()


gmx_check_if_changed(SIMD_CHANGED GMX_SIMD_ACTIVE)
if (SIMD_CHANGED AND DEFINED SIMD_STATUS_MESSAGE)
    message(STATUS "${SIMD_STATUS_MESSAGE}")
endif()

# By default, 32-bit windows cannot pass SIMD (SSE/AVX) arguments in registers,
# and even on 64-bit (all platforms) it is only used for a handful of arguments.
# The __vectorcall (MSVC, from MSVC2013) or __regcall (ICC) calling conventions
# enable this, which is critical to enable 32-bit SIMD and improves performance
# for 64-bit SIMD.
# Check if the compiler supports one of these, and in that case set gmx_simdcall
# to that string. If we do not have any such calling convention modifier, set it
# to an empty string.
#
# Update 2015-11-04: As of version 3.6, clang has added support for __vectorcall
# (also on Linux). This appears to be buggy for the reference SIMD
# implementation when using the Debug build (when functions are not inlined) 
# while it seems works fine for the actual SIMD implementations. This is likely
# because the reference build ends up passing lots of structures with arrays
# rather than actual vector data. For now we disable __vectorcall with clang
# when using the reference build.
# 
# xlc 13.1.5 does not seem recognize any attribute, and warns about invalid ones
# so we avoid searching for any.
#
if(NOT DEFINED GMX_SIMD_CALLING_CONVENTION)
    if(GMX_TARGET_BGQ)
        set(CALLCONV_LIST " ")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND GMX_SIMD_ACTIVE STREQUAL "REFERENCE")
        set(CALLCONV_LIST __regcall " ")
   elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")
        set(CALLCONV_LIST " ")
    else()
        set(CALLCONV_LIST __vectorcall __regcall " ")
    endif()
    foreach(callconv ${CALLCONV_LIST})
        set(callconv_compile_var "_callconv_${callconv}")
        # Some compilers warn about targets for which attributes are
        # ignored (e.g. clang on ARM), and in such cases we want this
        # check to lead to using no attribute in subsequent GROMACS
        # compilation, to avoid issuing the warning for lots of files.
        check_c_source_compiles("
#pragma GCC diagnostic error \"-Wignored-attributes\"
int ${callconv} f(int i) {return i;} int main(void) {return f(0);}
" ${callconv_compile_var})
        if(${callconv_compile_var})
            set(GMX_SIMD_CALLING_CONVENTION "${callconv}" CACHE INTERNAL "Calling convention for SIMD routines" FORCE)
            break()
        endif()
    endforeach()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # GCC bug 49001, 54412 on Windows (just warn, since it might be fixed in later versions)
    if((CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.0" OR CMAKE_SIZEOF_VOID_P EQUAL 8)
            AND (WIN32 OR CYGWIN)
            AND (GMX_SIMD_ACTIVE MATCHES "AVX") AND NOT (GMX_SIMD_ACTIVE STREQUAL "AVX_128_FMA"))
        message(WARNING "GCC on Windows (GCC older than 4.9 in 32-bit mode, or any version in 64-bit mode) with 256-bit AVX will probably crash. You might want to choose a different GMX_SIMD or a different compiler.")
    endif()
endif()

endmacro()

