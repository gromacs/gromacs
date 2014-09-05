#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# include avx test source, used if the AVX flags are set below
include(gmxTestAVXMaskload)
include(gmxFindFlagsForSource)


macro(gmx_use_clang_as_with_gnu_compilers_on_osx)
    # On OS X, we often want to use gcc instead of clang, since gcc supports
    # OpenMP. However, by default gcc uses the external system assembler, which
    # does not support AVX, so we need to tell the linker to use the clang
    # compilers assembler instead - and this has to happen before we detect AVX
    # flags.
    if(APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
        gmx_test_cflag(GNU_C_USE_CLANG_AS "-Wa,-q" SIMD_C_FLAGS)
    endif()
    if(APPLE AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        gmx_test_cxxflag(GNU_CXX_USE_CLANG_AS "-Wa,-q" SIMD_CXX_FLAGS)
    endif()
endmacro()


macro(gmx_test_simd)
#
# To improve backward compatibility on x86 SIMD architectures,
# we set the flags for all SIMD instructions that are supported, not only
# the most recent instruction set. I.e., if your machine supports AVX2_256,
# we will set flags both for AVX2_256, AVX_256, SSE4.1, and SSE2 support.

if(${GMX_SIMD} STREQUAL "NONE")
    # nothing to do configuration-wise
    set(SIMD_STATUS_MESSAGE "SIMD instructions disabled")
elseif(${GMX_SIMD} STREQUAL "SSE2")

    gmx_find_cflag_for_source(CFLAGS_SSE2 "C compiler SSE2 flag"
                              "#include<xmmintrin.h>
                              int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_rsqrt_ps(x);return _mm_movemask_ps(x);}"
                              SIMD_C_FLAGS
                              "-msse2" "/arch:SSE2" "-hgnu")
    gmx_find_cxxflag_for_source(CXXFLAGS_SSE2 "C++ compiler SSE2 flag"
                                "#include<xmmintrin.h>
                                int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_rsqrt_ps(x);return _mm_movemask_ps(x);}"
                                SIMD_CXX_FLAGS
                                "-msse2" "/arch:SSE2" "-hgnu")

    if(NOT CFLAGS_SSE2 OR NOT CXXFLAGS_SSE2)
        message(FATAL_ERROR "Cannot find SSE2 compiler flag. Use a newer compiler, or disable SIMD (slower).")
    endif()

    set(GMX_SIMD_X86_SSE2 1)
    set(SIMD_STATUS_MESSAGE "Enabling SSE2 SIMD instructions")

elseif(${GMX_SIMD} STREQUAL "SSE4.1")

    # Note: MSVC enables SSE4.1 with the SSE2 flag, so we include that in testing.
    gmx_find_cflag_for_source(CFLAGS_SSE4_1 "C compiler SSE4.1 flag"
                              "#include<smmintrin.h>
                              int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_dp_ps(x,x,0x77);return _mm_movemask_ps(x);}"
                              SIMD_C_FLAGS
                              "-msse4.1" "/arch:SSE4.1" "/arch:SSE2" "-hgnu")
    gmx_find_cxxflag_for_source(CXXFLAGS_SSE4_1 "C++ compiler SSE4.1 flag"
                                "#include<smmintrin.h>
                                int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_dp_ps(x,x,0x77);return _mm_movemask_ps(x);}"
                                SIMD_CXX_FLAGS
                                "-msse4.1" "/arch:SSE4.1" "/arch:SSE2" "-hgnu")

    if(NOT CFLAGS_SSE4_1 OR NOT CXXFLAGS_SSE4_1)
        message(FATAL_ERROR "Cannot find SSE4.1 compiler flag. "
                            "Use a newer compiler, or choose SSE2 SIMD (slower).")
    endif()

    if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND CMAKE_C_COMPILER_VERSION VERSION_EQUAL "11.1")
        message(FATAL_ERROR "You are using Intel compiler version 11.1, which produces incorrect results with SSE4.1 SIMD. You need to use a newer compiler (e.g. icc >= 12.0) or in worst case try a lower level of SIMD if performance is not critical.")
    endif()

    set(GMX_SIMD_X86_SSE4_1 1)
    set(SIMD_STATUS_MESSAGE "Enabling SSE4.1 SIMD instructions")

elseif(${GMX_SIMD} STREQUAL "AVX_128_FMA")

    gmx_use_clang_as_with_gnu_compilers_on_osx()

    # AVX128/FMA on AMD is a bit complicated. We need to do detection in three stages:
    # 1) Find the flags required for generic AVX support
    # 2) Find the flags necessary to enable fused-multiply add support
    # 3) Optional: Find a flag to enable the AMD XOP instructions

    ### STAGE 1: Find the generic AVX flag
    gmx_find_cflag_for_source(CFLAGS_AVX_128 "C compiler AVX (128 bit) flag"
                              "#include<immintrin.h>
                              int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_permute_ps(x,1);return 0;}"
                              SIMD_C_FLAGS
                              "-mavx" "/arch:AVX" "-hgnu")
    gmx_find_cxxflag_for_source(CXXFLAGS_AVX_128 "C++ compiler AVX (128 bit) flag"
                                "#include<immintrin.h>
                                int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_permute_ps(x,1);return 0;}"
                                SIMD_CXX_FLAGS
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

    gmx_find_cflag_for_source(CFLAGS_AVX_128_FMA "C compiler AVX (128 bit) FMA4 flag"
"#include<immintrin.h>
${INCLUDE_X86INTRIN_H}
${INCLUDE_INTRIN_H}
int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_macc_ps(x,x,x);return _mm_movemask_ps(x);}"
                              SIMD_C_FLAGS
                              "-mfma4" "-hgnu")
    gmx_find_cxxflag_for_source(CXXFLAGS_AVX_128_FMA "C++ compiler AVX (128 bit) FMA4 flag"
"#include<immintrin.h>
${INCLUDE_X86INTRIN_H}
${INCLUDE_INTRIN_H}
int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_macc_ps(x,x,x);return _mm_movemask_ps(x);}"
                                SIMD_CXX_FLAGS
                                "-mfma4" "-hgnu")

    # We only need to check the last (FMA) test; that will always fail if the basic AVX128 test failed
    if(NOT CFLAGS_AVX_128_FMA OR NOT CXXFLAGS_AVX_128_FMA)
        message(FATAL_ERROR "Cannot find compiler flags for 128 bit AVX with FMA support. Use a newer compiler, or choose SSE4.1 SIMD (slower).")
    endif()

    ### STAGE 3: Optional: Find the XOP instruction flag (No point in yelling if this does not work)
    gmx_find_cflag_for_source(CFLAGS_AVX_128_XOP "C compiler AVX (128 bit) XOP flag"
"#include<immintrin.h>
${INCLUDE_X86INTRIN_H}
${INCLUDE_INTRIN_H}
int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_frcz_ps(x);return _mm_movemask_ps(x);}"
                              SIMD_C_FLAGS
                              "-mxop")
    gmx_find_cxxflag_for_source(CXXFLAGS_AVX_128_XOP "C++ compiler AVX (128 bit) XOP flag"
"#include<immintrin.h>
${INCLUDE_X86INTRIN_H}
${INCLUDE_INTRIN_H}
int main(){__m128 x=_mm_set1_ps(0.5);x=_mm_frcz_ps(x);return _mm_movemask_ps(x);}"
                                SIMD_CXX_FLAGS
                                "-mxop")

    # We don't have the full compiler version string yet (BUILD_C_COMPILER),
    # so we can't distinguish vanilla from Apple clang versions, but catering for a few rare AMD
    # hackintoshes is not worth the effort.
    if (APPLE AND ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang" OR
                "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
        message(WARNING "Due to a known compiler bug, Clang up to version 3.2 (and Apple Clang up to version 4.1) produces incorrect code with AVX_128_FMA SIMD. As we cannot work around this bug on OS X, you will have to select a different compiler or SIMD instruction set.")
    endif()


    if (GMX_USE_CLANG_C_FMA_BUG_WORKAROUND)
        # we assume that we have an external assembler that supports AVX
        message(STATUS "Clang ${CMAKE_C_COMPILER_VERSION} detected, enabling FMA bug workaround")
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -no-integrated-as")
    endif()
    if (GMX_USE_CLANG_CXX_FMA_BUG_WORKAROUND)
        # we assume that we have an external assembler that supports AVX
        message(STATUS "Clang ${CMAKE_CXX_COMPILER_VERSION} detected, enabling FMA bug workaround")
        set(EXTRA_CXX_FLAGS "${EXTRA_CXX_FLAGS} -no-integrated-as")
    endif()

    gmx_test_avx_gcc_maskload_bug(GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG "${SIMD_C_FLAGS}")

    set(GMX_SIMD_X86_AVX_128_FMA 1)
    set(SIMD_STATUS_MESSAGE "Enabling 128-bit AVX SIMD Gromacs SIMD (with fused-multiply add)")

elseif(${GMX_SIMD} STREQUAL "AVX_256")

    gmx_use_clang_as_with_gnu_compilers_on_osx()

    gmx_find_cflag_for_source(CFLAGS_AVX "C compiler AVX (256 bit) flag"
                              "#include<immintrin.h>
                              int main(){__m256 x=_mm256_set1_ps(0.5);x=_mm256_add_ps(x,x);return _mm256_movemask_ps(x);}"
                              SIMD_C_FLAGS
                              "-mavx" "/arch:AVX" "-hgnu")
    gmx_find_cxxflag_for_source(CXXFLAGS_AVX "C++ compiler AVX (256 bit) flag"
                                "#include<immintrin.h>
                                int main(){__m256 x=_mm256_set1_ps(0.5);x=_mm256_add_ps(x,x);return _mm256_movemask_ps(x);}"
                                SIMD_CXX_FLAGS
                                "-mavx" "/arch:AVX" "-hgnu")

    if(NOT CFLAGS_AVX OR NOT CXXFLAGS_AVX)
        message(FATAL_ERROR "Cannot find AVX compiler flag. Use a newer compiler, or choose SSE4.1 SIMD (slower).")
    endif()

    gmx_test_avx_gcc_maskload_bug(GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG "${SIMD_C_FLAGS}")

    set(GMX_SIMD_X86_AVX_256 1)
    set(SIMD_STATUS_MESSAGE "Enabling 256-bit AVX SIMD instructions")

elseif(${GMX_SIMD} STREQUAL "AVX2_256")

    gmx_use_clang_as_with_gnu_compilers_on_osx()

    gmx_find_cflag_for_source(CFLAGS_AVX2 "C compiler AVX2 flag"
                              "#include<immintrin.h>
                              int main(){__m256i x=_mm256_set1_epi32(5);x=_mm256_add_epi32(x,x);return _mm256_movemask_epi8(x);}"
                              SIMD_C_FLAGS
                              "-march=core-avx2" "-mavx2" "/arch:AVX" "-hgnu") # no AVX2-specific flag for MSVC yet
    gmx_find_cxxflag_for_source(CXXFLAGS_AVX2 "C++ compiler AVX2 flag"
                                "#include<immintrin.h>
                                int main(){__m256i x=_mm256_set1_epi32(5);x=_mm256_add_epi32(x,x);return _mm256_movemask_epi8(x);}"
                                SIMD_CXX_FLAGS
                                "-march=core-avx2" "-mavx2" "/arch:AVX" "-hgnu") # no AVX2-specific flag for MSVC yet

    if(NOT CFLAGS_AVX2 OR NOT CXXFLAGS_AVX2)
        message(FATAL_ERROR "Cannot find AVX2 compiler flag. Use a newer compiler, or choose AVX SIMD (slower).")
    endif()

    # No need to test for Maskload bug - it was fixed before gcc added AVX2 support

    set(GMX_SIMD_X86_AVX2_256 1)
    set(SIMD_STATUS_MESSAGE "Enabling 256-bit AVX2 SIMD instructions")

elseif(${GMX_SIMD} STREQUAL "IBM_QPX")

    try_compile(TEST_QPX ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestQPX.c")

    if (TEST_QPX)
        message(WARNING "IBM QPX SIMD instructions selected. This will work, but SIMD kernels are only available for the Verlet cut-off scheme. The plain C kernels that are used for the group cut-off scheme kernels will be slow, so please consider using the Verlet cut-off scheme.")
        set(GMX_SIMD_IBM_QPX 1)
        set(SIMD_STATUS_MESSAGE "Enabling IBM QPX SIMD instructions")

    else()
        message(FATAL_ERROR "Cannot compile the requested IBM QPX intrinsics. If you are compiling for BlueGene/Q with the XL compilers, use 'cmake .. -DCMAKE_TOOLCHAIN_FILE=Platform/BlueGeneQ-static-XL-C' to set up the tool chain.")
    endif()

elseif(${GMX_SIMD} STREQUAL "SPARC64_HPC_ACE")

    # Note that GMX_RELAXED_DOUBLE_PRECISION is enabled by default in the top-level CMakeLists.txt

    set(GMX_SIMD_SPARC64_HPC_ACE 1)
    set(SIMD_STATUS_MESSAGE "Enabling Sparc64 HPC-ACE SIMD instructions")

elseif(${GMX_SIMD} STREQUAL "REFERENCE")

    # NB: This file handles settings for the SIMD module, so in the interest 
    # of proper modularization, please do NOT put any verlet kernel settings in this file.

    if(GMX_SIMD_REF_FLOAT_WIDTH)
        add_definitions(-DGMX_SIMD_REF_FLOAT_WIDTH=${GMX_SIMD_REF_FLOAT_WIDTH})
    endif()
    if(GMX_SIMD_REF_DOUBLE_WIDTH)
      	add_definitions(-DGMX_SIMD_REF_DOUBLE_WIDTH=${GMX_SIMD_REF_DOUBLE_WIDTH})
    endif()

    set(GMX_SIMD_REFERENCE 1)
    set(SIMD_STATUS_MESSAGE "Enabling reference (emulated) SIMD instructions.")

else()
    gmx_invalid_option_value(GMX_SIMD)
endif()


gmx_check_if_changed(SIMD_CHANGED GMX_SIMD)
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
if(NOT DEFINED GMX_SIMD_CALLING_CONVENTION)
    foreach(callconv __vectorcall __regcall "")
        set(callconv_compile_var "_callconv_${callconv}")
        check_c_source_compiles("int ${callconv} f(int i) {return i;} int main(void) {return f(0);}" ${callconv_compile_var})
        if(${callconv_compile_var})
            set(GMX_SIMD_CALLING_CONVENTION "${callconv}" CACHE INTERNAL "Calling convention for SIMD routines" FORCE)
            break()
        endif()
    endforeach()
endif()

endmacro()

