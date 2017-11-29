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
include(gmxSimdFlags)

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

    gmx_find_simd_sse2_flags(SIMD_SSE2_C_SUPPORTED SIMD_SSE2_CXX_SUPPORTED
                             SIMD_SSE2_C_FLAGS SIMD_SSE2_CXX_FLAGS)

    if(NOT SIMD_SSE2_C_SUPPORTED OR NOT SIMD_SSE2_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("SSE2" "disable SIMD support (slow)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_SSE2_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_SSE2_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling SSE2 SIMD instructions using CXX flags: ${SIMD_SSE2_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "SSE4.1")

    gmx_find_simd_sse4_1_flags(SIMD_SSE4_1_C_SUPPORTED SIMD_SSE4_1_CXX_SUPPORTED
                               SIMD_SSE4_1_C_FLAGS SIMD_SSE4_1_CXX_FLAGS)

    if(NOT SIMD_SSE4_1_C_SUPPORTED OR NOT SIMD_SSE4_1_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("SSE4.1" "choose SSE2 SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_SSE4_1_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_SSE4_1_CXX_FLAGS}")
    set(GMX_SIMD_X86_SSE4_1 1)
    set(SIMD_STATUS_MESSAGE "Enabling SSE4.1 SIMD instructions using CXX flags: ${SIMD_SSE4_1_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_128_FMA")

    gmx_find_simd_avx_128_fma_flags(SIMD_AVX_128_FMA_C_SUPPORTED SIMD_AVX_128_FMA_CXX_SUPPORTED
                                    SIMD_AVX_128_FMA_C_FLAGS SIMD_AVX_128_FMA_CXX_FLAGS)

    if(NOT SIMD_AVX_128_FMA_C_SUPPORTED OR NOT SIMD_AVX_128_FMA_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("128-bit AVX with FMA support" "choose SSE4.1 SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_AVX_128_FMA_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_AVX_128_FMA_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 128-bit AMD FMA SIMD instructions using CXX flags: ${SIMD_AVX_128_FMA_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_256")

    gmx_find_simd_avx_flags(SIMD_AVX_C_SUPPORTED SIMD_AVX_CXX_SUPPORTED
                            SIMD_AVX_C_FLAGS SIMD_AVX_CXX_FLAGS)

    if(NOT SIMD_AVX_C_SUPPORTED OR NOT SIMD_AVX_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("AVX" "choose SSE4.1 SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_AVX_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_AVX_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 256-bit AVX SIMD instructions using CXX flags: ${SIMD_AVX_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE MATCHES "AVX2_")

    gmx_find_simd_avx2_flags(SIMD_AVX2_C_SUPPORTED SIMD_AVX2_CXX_SUPPORTED
                             SIMD_AVX2_C_FLAGS SIMD_AVX2_CXX_FLAGS)

    if(NOT SIMD_AVX2_C_SUPPORTED OR NOT SIMD_AVX2_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("AVX2" "choose AVX SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_AVX2_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_AVX2_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)

    if(GMX_SIMD_ACTIVE STREQUAL "AVX2_128")
        set(SIMD_STATUS_MESSAGE "Enabling 128-bit AVX2 SIMD instructions using CXX flags: ${SIMD_AVX2_CXX_FLAGS}")
    else()
        set(SIMD_STATUS_MESSAGE "Enabling 256-bit AVX2 SIMD instructions using CXX flags: ${SIMD_AVX2_CXX_FLAGS}")
    endif()

elseif(GMX_SIMD_ACTIVE STREQUAL "MIC")

    # No flags needed. Not testing.
    set(GMX_SIMD_X86_MIC 1)
    set(SIMD_STATUS_MESSAGE "Enabling MIC (Xeon Phi) SIMD instructions without special flags.")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_512")

    gmx_find_simd_avx_512_flags(SIMD_AVX_512_C_SUPPORTED SIMD_AVX_512_CXX_SUPPORTED
                                SIMD_AVX_512_C_FLAGS SIMD_AVX_512_CXX_FLAGS)

    if(NOT SIMD_AVX_512_C_SUPPORTED OR NOT SIMD_AVX_512_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("AVX 512F" "choose a lower level of SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_AVX_512_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_AVX_512_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 512-bit AVX-512 SIMD instructions using CXX flags: ${SIMD_AVX_512_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "AVX_512_KNL")

    gmx_find_simd_avx_512_knl_flags(SIMD_AVX_512_KNL_C_SUPPORTED SIMD_AVX_512_KNL_CXX_SUPPORTED
                                    SIMD_AVX_512_KNL_C_FLAGS SIMD_AVX_512_KNL_CXX_FLAGS)

    if(NOT SIMD_AVX_512_KNL_C_SUPPORTED OR NOT SIMD_AVX_512_KNL_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("AVX 512ER" "choose a lower level of SIMD (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_AVX_512_KNL_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_AVX_512_KNL_CXX_FLAGS}")
    set(GMX_SIMD_X86_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 512-bit AVX-512-KNL SIMD instructions using CXX flags: ${SIMD_AVX_512_KNL_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "ARM_NEON")

    if (GMX_DOUBLE)
        message(FATAL_ERROR "ARM_NEON SIMD support is not available for a double precision build because the architecture lacks double-precision support")
    endif()

    gmx_find_simd_arm_neon_flags(SIMD_ARM_NEON_C_SUPPORTED SIMD_ARM_NEON_CXX_SUPPORTED
                                 SIMD_ARM_NEON_C_FLAGS SIMD_ARM_NEON_CXX_FLAGS)

    if(NOT SIMD_ARM_NEON_C_SUPPORTED OR NOT SIMD_ARM_NEON_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("ARM NEON" "disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_ARM_NEON_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_ARM_NEON_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling 32-bit ARM NEON SIMD instructions using CXX flags: ${SIMD_ARM_NEON_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "ARM_NEON_ASIMD")

    gmx_find_simd_arm_neon_asimd_flags(SIMD_ARM_NEON_ASIMD_C_SUPPORTED SIMD_ARM_NEON_ASIMD_CXX_SUPPORTED
                                       SIMD_ARM_NEON_ASIMD_C_FLAGS SIMD_ARM_NEON_ASIMD_CXX_FLAGS)

    if(NOT SIMD_ARM_NEON_ASIMD_C_SUPPORTED OR NOT SIMD_ARM_NEON_ASIMD_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("ARM (AArch64) NEON Advanced SIMD" "particularly gcc version 4.9 or later, or disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_ARM_NEON_ASIMD_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_ARM_NEON_ASIMD_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling ARM (AArch64) NEON Advanced SIMD instructions using CXX flags: ${SIMD_ARM_NEON_ASIMD_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "IBM_QPX")

    try_compile(TEST_QPX ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestQPX.c")

    if (TEST_QPX)
        set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
        set(SIMD_STATUS_MESSAGE "Enabling IBM QPX SIMD instructions without special flags.")
    else()
        gmx_give_fatal_error_when_simd_support_not_found("IBM QPX" "or 'cmake .. -DCMAKE_TOOLCHAIN_FILE=Platform/BlueGeneQ-static-bgclang-CXX' to set up the tool chain" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

elseif(GMX_SIMD_ACTIVE STREQUAL "IBM_VMX")

    gmx_find_simd_ibm_vmx_flags(SIMD_IBM_VMX_C_SUPPORTED SIMD_IBM_VMX_CXX_SUPPORTED
                                SIMD_IBM_VMX_C_FLAGS SIMD_IBM_VMX_CXX_FLAGS)

    if(NOT SIMD_IBM_VMX_C_SUPPORTED OR NOT SIMD_IBM_VMX_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("IBM VMX" "disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_IBM_VMX_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_IBM_VMX_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling IBM VMX SIMD instructions using CXX flags: ${SIMD_IBM_VMX_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "IBM_VSX")

    gmx_find_simd_ibm_vsx_flags(SIMD_IBM_VSX_C_SUPPORTED SIMD_IBM_VSX_CXX_SUPPORTED
                                SIMD_IBM_VSX_C_FLAGS SIMD_IBM_VSX_CXX_FLAGS)

    # Usually we check also for the C compiler here, but a C compiler
    # is not required for SIMD support on this platform. cmake through
    # at least version 3.7 cannot pass this check with the C compiler
    # in the latest xlc 13.1.5, but the C++ compiler has different
    # behaviour and is OK. See Redmine #2102.
    if(NOT SIMD_IBM_VSX_CXX_SUPPORTED)
        gmx_give_fatal_error_when_simd_support_not_found("IBM VSX" "disable SIMD support (slower)" "${SUGGEST_BINUTILS_UPDATE}")
    endif()

    set(SIMD_C_FLAGS "${SIMD_IBM_VSX_C_FLAGS}")
    set(SIMD_CXX_FLAGS "${SIMD_IBM_VSX_CXX_FLAGS}")
    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling IBM VSX SIMD instructions using CXX flags: ${SIMD_IBM_VSX_CXX_FLAGS}")

elseif(GMX_SIMD_ACTIVE STREQUAL "SPARC64_HPC_ACE")

    # Note that GMX_RELAXED_DOUBLE_PRECISION is enabled by default in the top-level CMakeLists.txt

    set(GMX_SIMD_${GMX_SIMD_ACTIVE} 1)
    set(SIMD_STATUS_MESSAGE "Enabling Sparc64 HPC-ACE SIMD instructions without special flags.")

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
    set(SIMD_STATUS_MESSAGE "Enabling reference (emulated) SIMD instructions without special flags.")

else()
    gmx_invalid_option_value(GMX_SIMD_ACTIVE)
endif()


gmx_check_if_changed(SIMD_CHANGED GMX_SIMD_ACTIVE)
if (SIMD_CHANGED AND DEFINED SIMD_STATUS_MESSAGE)
    message(STATUS "${SIMD_STATUS_MESSAGE}")
endif()

# While AVX-512 is a more recent SIMD ISA than AVX2, some Intel CPUs only have
# a single AVX-512 FMA unit, but two AVX2 FMA units, and then it is better to
# use AVX2. The only way to test this is to execute a small timing loop.
# To be able to recommend the user whether s/he should try AVX-512 instead of
# AVX2, we need to compile a single file with AVX512 flags. We do this
# automatically, but this option provides a way to turn it off in case it
# breaks something. The actual test source file is built if
# SIMD_AVX_512_CXX_SUPPORTED is set, so it will always be included if we have
# GMX_SIMD=AVX_512.
set(GMX_ENABLE_AVX512_TESTS ON CACHE INTERNAL "Compile AVX512 code to test FMA units, even when not using AVX512 SIMD")
mark_as_advanced(GMX_ENABLE_AVX512_TESTS)

if(GMX_ENABLE_AVX512_TESTS AND
    (GMX_SIMD_ACTIVE STREQUAL "AVX_256" OR GMX_SIMD_ACTIVE STREQUAL "AVX2_256" OR GMX_SIMD_ACTIVE STREQUAL "AVX2_128"))
    if(NOT DEFINED SIMD_AVX_512_CXX_SUPPORTED)
        message(STATUS "Detecting flags to enable runtime detection of AVX-512 units on newer CPUs")
        set(SIMD_AVX_512_REPORT_STATUS 1)
    endif()
    gmx_find_simd_avx_512_flags(SIMD_AVX_512_C_SUPPORTED SIMD_AVX_512_CXX_SUPPORTED
                                SIMD_AVX_512_C_FLAGS SIMD_AVX_512_CXX_FLAGS)
    if(SIMD_AVX_512_REPORT_STATUS)
        if(SIMD_AVX_512_CXX_SUPPORTED)
            message(STATUS "Detecting flags to enable runtime detection of AVX-512 units on newer CPUs - ${SIMD_AVX_512_CXX_FLAGS}")
        else()
            message(STATUS "Detecting flags to enable runtime detection of AVX-512 units on newer CPUs - not supported")
        endif()
    endif()
    # Since we might be overriding AVX2 architecture flags with the AVX512 flags for the
    # files where it is used, we also check for a flag not to warn about the first (unused) arch.
    # To avoid spamming the user with lots of gromacs tests we just call the CMake flag test directly.
    foreach(_testflag "-Wno-unused-command-line-argument" "-wd10121")
        string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" FLAG_ACCEPTED_VARIABLE "${_testflag}_FLAG_ACCEPTED")
        check_cxx_compiler_flag("${_testflag}" ${FLAG_ACCEPTED_VARIABLE})
        if(${FLAG_ACCEPTED_VARIABLE})
            set(CXX_NO_UNUSED_OPTION_WARNING_FLAGS "${_testflag}")
            break()
        endif()
    endforeach(_testflag)
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

