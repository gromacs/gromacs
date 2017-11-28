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

include(gmxFindFlagsForSource)

# Macro that manages setting the respective C and C++ toolchain
# variables so that subsequent tests for SIMD support can work.
macro(find_x86_toolchain_flags TOOLCHAIN_C_FLAGS_VARIABLE TOOLCHAIN_CXX_FLAGS_VARIABLE)
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
macro(find_power_vsx_toolchain_flags TOOLCHAIN_C_FLAGS_VARIABLE TOOLCHAIN_CXX_FLAGS_VARIABLE)
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

# SSE2
function(gmx_find_simd_sse2_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)
    gmx_find_flags(SIMD_SSE2_C_FLAGS_RESULT SIMD_SSE2_CXX_FLAGS_RESULT
        "#include<xmmintrin.h>
         int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_rsqrt_ps(x);return _mm_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_SSE2_C_FLAGS SIMD_SSE2_CXX_FLAGS
        "-msse2" "/arch:SSE2" "-hgnu")

    if(${SIMD_SSE2_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_SSE2_C_FLAGS}" CACHE INTERNAL "C flags required for SSE2 instructions")
    endif()
    if(${SIMD_SSE2_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_SSE2_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for SSE2 instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_SSE2_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for SSE2 C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_SSE2_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for SSE2 C++ flags" FORCE)
endfunction()

# SSE4.1
function(gmx_find_simd_sse4_1_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)
    # Note: MSVC enables SSE4.1 with the SSE2 flag, so we include that in testing.
    gmx_find_flags(SIMD_SSE4_1_C_FLAGS_RESULT SIMD_SSE4_1_CXX_FLAGS_RESULT
        "#include<smmintrin.h>
        int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_dp_ps(x,x,0x77);return _mm_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_SSE4_1_C_FLAGS SIMD_SSE4_1_CXX_FLAGS
        "-msse4.1" "/arch:SSE4.1" "/arch:SSE2" "-hgnu")

    if(${SIMD_SSE4_1_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_SSE4_1_C_FLAGS}" CACHE INTERNAL "C flags required for SSE4.1 instructions")
    endif()
    if(${SIMD_SSE4_1_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_SSE4_1_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for SSE4.1 instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_SSE4_1_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for SSE4.1 C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_SSE4_1_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for SSE4.1 C++ flags" FORCE)
endfunction()

# AVX, but using only 128-bit instructions and FMA (AMD XOP processors)
function(gmx_find_simd_avx_128_fma_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    # We don't have the full compiler version string yet (BUILD_C_COMPILER),
    # so we can't distinguish vanilla from Apple clang versions, but catering for a few rare AMD
    # hackintoshes is not worth the effort.
    if (APPLE AND (CMAKE_C_COMPILER_ID STREQUAL "Clang" OR
                CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
        message(WARNING "Due to a known compiler bug, Clang up to version 3.2 (and Apple Clang up to version 4.1) produces incorrect code with AVX_128_FMA SIMD.")
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

    ### STAGE 1: Find the generic AVX flag, but stick to 128-bit instructions
    gmx_find_flags(SIMD_AVX_128_FMA_C_FLAGS_RESULT SIMD_AVX_128_FMA_CXX_FLAGS_RESULT
        "#include<immintrin.h>
        int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_permute_ps(x,1);return 0;}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_AVX_GENERIC_C_FLAGS SIMD_AVX_GENERIC_CXX_FLAGS
        "-mavx" "/arch:AVX" "-hgnu")

    if(SIMD_AVX_128_FMA_C_FLAGS_RESULT AND SIMD_AVX_128_FMA_CXX_FLAGS_RESULT)
        set(MERGED_C_FLAGS "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX_GENERIC_C_FLAGS}")
        set(MERGED_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX_GENERIC_CXX_FLAGS}")

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

        gmx_find_flags(SIMD_AVX_128_FMA_C_FLAGS_RESULT SIMD_AVX_128_FMA_CXX_FLAGS_RESULT
            "#include<immintrin.h>
            ${INCLUDE_X86INTRIN_H}
            ${INCLUDE_INTRIN_H}
            int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_macc_ps(x,x,x);return _mm_movemask_ps(x);}"
            MERGED_C_FLAGS MERGED_CXX_FLAGS
            SIMD_AVX_AMD_FMA_C_FLAGS SIMD_AVX_AMD_FMA_CXX_FLAGS
            "-mfma4" "-hgnu")

        if(SIMD_AVX_128_FMA_C_FLAGS_RESULT AND SIMD_AVX_128_FMA_CXX_FLAGS_RESULT)
            set(MERGED_C_FLAGS "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX_AMD_FMA_C_FLAGS}")
            set(MERGED_CXX_FLAGS "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX_AMD_FMA_CXX_FLAGS}")
            ### STAGE 3: Find the XOP instruction flag. This is optional.
            gmx_find_flags(SIMD_AVX_XOP_C_FLAGS_RESULT SIMD_AVX_XOP_CXX_FLAGS_RESULT
                "#include<immintrin.h>
                ${INCLUDE_X86INTRIN_H}
                ${INCLUDE_INTRIN_H}
                int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_frcz_ps(x);return _mm_movemask_ps(x);}"
                MERGED_C_FLAGS MERGED_CXX_FLAGS
                SIMD_AVX_XOP_C_FLAGS SIMD_AVX_XOP_CXX_FLAGS
                "-mxop")
        endif()
    endif()

    if(${SIMD_AVX_128_FMA_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX_GENERIC_C_FLAGS} ${SIMD_AVX_AMD_FMA_C_FLAGS} ${SIMD_AVX_XOP_C_FLAGS}" CACHE INTERNAL "C flags required for 128-bit AVX with AMD FMA instructions")
    endif()
    if(${SIMD_AVX_128_FMA_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX_GENERIC_CXX_FLAGS} ${SIMD_AVX_AMD_FMA_CXX_FLAGS} ${SIMD_AVX_XOP_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for 128-bit AVX with AMD FMA instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_AVX_128_FMA_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for 128-bit AVX with AMD FMA C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_AVX_128_FMA_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for 128-bit AVX with AMD FMA C++ flags" FORCE)
endfunction()


# AVX (no AMD extensions)
function(gmx_find_simd_avx_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)
    gmx_find_flags(SIMD_AVX_C_FLAGS_RESULT SIMD_AVX_CXX_FLAGS_RESULT
        "#include<immintrin.h>
         int main(){__m256 x=_mm256_set1_ps(0.5);x=_mm256_add_ps(x,x);return _mm256_movemask_ps(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_AVX_C_FLAGS SIMD_AVX_CXX_FLAGS
        "-mavx" "/arch:AVX" "-hgnu")

    if(${SIMD_AVX_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX_C_FLAGS}" CACHE INTERNAL "C flags required for AVX instructions")
    endif()
    if(${SIMD_AVX_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for AVX instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_AVX_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_AVX_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX C++ flags" FORCE)
endfunction()

# AVX2
function(gmx_find_simd_avx2_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)
    gmx_find_flags(SIMD_AVX2_C_FLAGS_RESULT SIMD_AVX2_CXX_FLAGS_RESULT
        "#include<immintrin.h>
        int main(){__m256i x=_mm256_set1_epi32(5);x=_mm256_add_epi32(x,x);return _mm256_movemask_epi8(x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_AVX2_C_FLAGS SIMD_AVX2_CXX_FLAGS
        "-march=core-avx2" "-mavx2" "/arch:AVX" "-hgnu") # no AVX2-specific flag for MSVC yet

    if(${SIMD_AVX2_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX2_C_FLAGS}" CACHE INTERNAL "C flags required for AVX2 instructions")
    endif()
    if(${SIMD_AVX2_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX2_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for AVX2 instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_AVX2_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX2 C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_AVX2_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX2 C++ flags" FORCE)
endfunction()


# AVX-512F (Skylake-X)
function(gmx_find_simd_avx_512_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(SIMD_AVX_512_C_FLAGS_RESULT SIMD_AVX_512_CXX_FLAGS_RESULT
        "#include<immintrin.h>
         int main(){__m512 x=_mm512_set1_ps(0.5); __m512 y=_mm512_fmadd_ps(x,x,x);
          __m512i i = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
          __mmask16 mask = (short)(0xffff);
          int idata[16]; i  = _mm512_maskz_permutexvar_epi32(mask, i, i);
          _mm512_storeu_si512(idata, i);
          return idata[0]*(int)(_mm512_cmp_ps_mask(x,y,_CMP_LT_OS));}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_AVX_512_C_FLAGS SIMD_AVX_512_CXX_FLAGS
        "-xCORE-AVX512 -qopt-zmm-usage=high" "-xCORE-AVX512" "-mavx512f -mfma" "-mavx512f" "/arch:AVX" "-hgnu") # no AVX_512F flags known for MSVC yet. ICC should use ZMM if code anyhow uses ZMM

    if(${SIMD_AVX_512_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX_512_C_FLAGS}" CACHE INTERNAL "C flags required for AVX-512 instructions")
    endif()
    if(${SIMD_AVX_512_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX_512_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for AVX-512 instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_AVX_512_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX-512 C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_AVX_512_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX-512 C++ flags" FORCE)
endfunction()


# AVX-512ER (KNL)
function(gmx_find_simd_avx_512_knl_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_x86_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)

    gmx_find_flags(SIMD_AVX_512_KNL_C_FLAGS_RESULT SIMD_AVX_512_KNL_CXX_FLAGS_RESULT
        "#include<immintrin.h>
        int main(){__m512 y,x=_mm512_set1_ps(0.5);y=_mm512_rsqrt28_ps(x);return (int)_mm512_cmp_ps_mask(x,y,_CMP_LT_OS);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_AVX_512_KNL_C_FLAGS SIMD_AVX_512_KNL_CXX_FLAGS
        "-xMIC-AVX512" "-mavx512er -mfma" "-mavx512er" "/arch:AVX" "-hgnu") # no AVX_512ER flags known for MSVC yet

    if(${SIMD_AVX_512_KNL_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_AVX_512_KNL_C_FLAGS}" CACHE INTERNAL "C flags required for AVX-512 for KNL instructions")
    endif()
    if(${SIMD_AVX_512_KNL_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_AVX_512_KNL_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for AVX-512 for KNL instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_AVX_512_KNL_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX-512 for KNL C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_AVX_512_KNL_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for AVX-512 for KNL C++ flags" FORCE)
endfunction()


# Arm Neon (32-bit ARM)
function(gmx_find_simd_arm_neon_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)

    gmx_find_flags(SIMD_ARM_NEON_C_FLAGS_RESULT SIMD_ARM_NEON_CXX_FLAGS_RESULT
        "#include<arm_neon.h>
         int main(){float32x4_t x=vdupq_n_f32(0.5);x=vmlaq_f32(x,x,x);return vgetq_lane_f32(x,0)>0;}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_ARM_NEON_C_FLAGS SIMD_ARM_NEON_CXX_FLAGS
        "-mfpu=neon-vfpv4" "-mfpu=neon" "")

    if(${SIMD_ARM_NEON_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_ARM_NEON_C_FLAGS}" CACHE INTERNAL "C flags required for Arm Neon instructions")
    endif()
    if(${SIMD_ARM_NEON_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_ARM_NEON_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for Arm Neon instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_ARM_NEON_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for Arm Neon C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_ARM_NEON_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for Arm Neon C++ flags" FORCE)
endfunction()

# Arm Neon Asimd (64-bit ARM)
function(gmx_find_simd_arm_neon_asimd_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)

    gmx_find_flags(SIMD_ARM_NEON_ASIMD_C_FLAGS_RESULT SIMD_ARM_NEON_ASIMD_CXX_FLAGS_RESULT
        "#include<arm_neon.h>
         int main(){float64x2_t x=vdupq_n_f64(0.5);x=vfmaq_f64(x,x,x);x=vrndnq_f64(x);return vgetq_lane_f64(x,0)>0;}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_ARM_NEON_ASIMD_C_FLAGS SIMD_ARM_NEON_ASIMD_CXX_FLAGS
        "")

    if(${SIMD_ARM_NEON_ASIMD_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_ARM_NEON_ASIMD_C_FLAGS}" CACHE INTERNAL "C flags required for Arm Neon Asimd instructions")
    endif()
    if(${SIMD_ARM_NEON_ASIMD_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_ARM_NEON_ASIMD_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for Arm Neon Asimd instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_ARM_NEON_ASIMD_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for Arm Neon Asimd C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_ARM_NEON_ASIMD_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for Arm Neon Asimd C++ flags" FORCE)
endfunction()

# IBM VMX (power6)
function(gmx_find_simd_ibm_vmx_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)

    gmx_find_flags(SIMD_IBM_VMX_C_FLAGS_RESULT SIMD_IBM_VMX_CXX_FLAGS_RESULT
        "#include<altivec.h>
         int main(){vector float x,y=vec_ctf(vec_splat_s32(1),0);x=vec_madd(y,y,y);return vec_all_ge(y,x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_IBM_VMX_C_FLAGS SIMD_IBM_VMX_CXX_FLAGS
        "-maltivec -mabi=altivec" "-qarch=auto -qaltivec")

    if(${SIMD_IBM_VMX_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_IBM_VMX_C_FLAGS}" CACHE INTERNAL "C flags required for IBM VMX instructions")
    endif()
    if(${SIMD_IBM_VMX_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_IBM_VMX_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for IBM VMX instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_IBM_VMX_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for IBM VMX C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_IBM_VMX_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for IBM VMX C++ flags" FORCE)
endfunction()

# IBM VSX (power7 and later)
function(gmx_find_simd_ibm_vsx_flags C_FLAGS_RESULT CXX_FLAGS_RESULT C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    find_power_vsx_toolchain_flags(TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS)
    gmx_find_flags(SIMD_IBM_VSX_C_FLAGS_RESULT SIMD_IBM_VSX_CXX_FLAGS_RESULT
        "#include<altivec.h>
         int main(){vector double x,y=vec_splats(1.0);x=vec_madd(y,y,y);return vec_all_ge(y,x);}"
        TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS
        SIMD_IBM_VSX_C_FLAGS SIMD_IBM_VSX_CXX_FLAGS
        "-mvsx" "-maltivec -mabi=altivec" "-qarch=auto -qaltivec")

    if(${SIMD_IBM_VSX_C_FLAGS_RESULT})
        set(${C_FLAGS_VARIABLE} "${TOOLCHAIN_C_FLAGS} ${SIMD_IBM_VSX_C_FLAGS}" CACHE INTERNAL "C flags required for IBM VSX instructions")
    endif()
    if(${SIMD_IBM_VSX_CXX_FLAGS_RESULT})
        set(${CXX_FLAGS_VARIABLE} "${TOOLCHAIN_CXX_FLAGS} ${SIMD_IBM_VSX_CXX_FLAGS}" CACHE INTERNAL "C++ flags required for IBM VSX instructions")
    endif()
    set(${C_FLAGS_RESULT} ${SIMD_IBM_VSX_C_FLAGS_RESULT} CACHE INTERNAL "Result of test for IBM VSX C flags" FORCE)
    set(${CXX_FLAGS_RESULT} ${SIMD_IBM_VSX_CXX_FLAGS_RESULT} CACHE INTERNAL "Result of test for IBM VSX C++ flags" FORCE)
endfunction()
