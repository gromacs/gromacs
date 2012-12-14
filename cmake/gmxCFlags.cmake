#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#

# Test C flags FLAGS, and set cache VARIABLE for the flags if they
# work. Also add the flags to CFLAGSVAR (which should not be
# cached). Previous values in VARIABLE (from earlier iterations or the
# user) are preserved.
macro(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)
    # Run the test only if we don't know whether the flags might work.
    if(NOT DEFINED ${VARIABLE})
        check_c_compiler_flag("${FLAGS}" _testresult)
        if(_testresult)
            # If the test was conducted and successful, preserve the
            # result of the test in the cache.
            SET(${VARIABLE} ${FLAGS} CACHE STRING "" FORCE)
        endif()
    endif()

    # At this point, if ${VARIABLE} exists, it contains flags that
    # CMake has tested or the user has provided, and which will stay
    # in the cache. (Though if the user didn't specify a variable type
    # on the command line, they might not see it in ccmake.)
    if(${VARIABLE})
        # Users do not need to know about test results
        mark_as_advanced(${VARIABLE})
        # Add the flags to an (uncached) list of flags to be used in
        # this iteration of CMake
        set(${CFLAGSVAR} "${${VARIABLE}} ${${CFLAGSVAR}}")
    endif(${VARIABLE})
    # If ${VARIABLE} does not exist, the calling routine should handle
    # communicating back to the user.
endmacro()

# Test C++ flags FLAGS, and set cache VARIABLE for the flags if they
# work. Also add the flags to CXXFLAGSVAR (which should not be
# cached). Previous values in VARIABLE (from earlier iterations or the
# user) are preserved.
macro(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)
    # Run the test only if we don't know whether the flags might work.
    if(NOT DEFINED ${VARIABLE})
        check_cxx_compiler_flag("${FLAGS}" _testresult)
        if(_testresult)
            # If the test was conducted and successful, preserve the
            # result of the test in the cache.
            SET(${VARIABLE} ${FLAGS} CACHE STRING "" FORCE)
        endif()
    endif()

    # At this point, if ${VARIABLE} exists, it contains flags that
    # CMake has tested or the user has provided, and which will stay
    # in the cache. (Though if the user didn't specify a variable type
    # on the command line, they might not see it in ccmake.)
    if(${VARIABLE})
        # Users do not need to know about test results
        mark_as_advanced(${VARIABLE})
        # Add the flags to an (uncached) list of flags to be used in
        # this iteration of CMake
        set(${CXXFLAGSVAR} "${${VARIABLE}} ${${CXXFLAGSVAR}}")
    endif(${VARIABLE})
    # If ${VARIABLE} does not exist, the calling routine should handle
    # communicating back to the user.
endmacro()


# This is the actual exported function to be called 
MACRO(gmx_c_flags)

    include(CheckCCompilerFlag)
    include(CheckCXXCompilerFlag)

    # First, set up general compiler optimization flags for different
    # compilers and build types (Release and Debug)

    # Note that caching the results of tests is worthwhile, because
    # the only dependency they should have is the compiler, and that
    # cannot change. Caching combinations of test results leads to
    # problems with consistency, currency and brevity, so we rebuild
    # that on the fly as needed.

    # gcc
    if(CMAKE_COMPILER_IS_GNUCC)

        #Fix for LLVM OpenMP bug (redmine 900). Needs to run before OpenMP flags are set below.
        if(GMX_OPENMP)
            exec_program(${CMAKE_C_COMPILER} ARGS --version OUTPUT_VARIABLE _compiler_output)
            if(_compiler_output MATCHES "llvm.*4\\.2")
                message(STATUS "OpenMP multithreading not supported with llvm-gcc 4.2, disabled")
                set(GMX_OPENMP OFF CACHE BOOL
                    "OpenMP multithreading not not supported with llvm-gcc 4.2, disabled!" FORCE)
            endif()
        endif()

        #flags are added in reverse order and -Wno* need to appear after -Wall
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare" GMXC_CFLAGS)
        # new in gcc 4.5
        GMX_TEST_CFLAG(CFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_COPT "-fomit-frame-pointer -funroll-all-loops"
                       GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_NOINLINE "-fno-inline" GMXC_CFLAGS_DEBUG)
    endif()
    # g++
    if(CMAKE_COMPILER_IS_GNUCXX)
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare" GMXC_CXXFLAGS)
      # new in gcc 4.5
        GMX_TEST_CXXFLAG(CXXFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_COPT "-fomit-frame-pointer -funroll-all-loops"
                         GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_NOINLINE "-fno-inline" GMXC_CXXFLAGS_DEBUG)
    endif()

    # icc
    if (CMAKE_C_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_STDGNU "-std=gnu99" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_OPT "-ip -funroll-all-loops" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_SSE2 "-msse2" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_X86 "-mtune=core2" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_IA64 "-mtune=itanium2" GMXC_CFLAGS_RELEASE)
        else()
            GMX_TEST_CFLAG(CFLAGS_WARN "/W2" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_SSE2 "/arch:SSE2" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_X86 "/Qip" GMXC_CFLAGS_RELEASE)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-ip -funroll-all-loops" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_SSE2 "-msse2" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_X86 "-mtune=core2" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_IA64 "-mtune=itanium2" 
                              GMXC_CXXFLAGS_RELEASE)
        else()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/W2" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_SSE2 "/arch:SSE2" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_X86 "/Qip" GMXC_CXXFLAGS_RELEASE)
        endif()
    endif()

    # pgi
    if (CMAKE_C_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CFLAG(CFLAGS_OPT "-fastsse" GMXC_CFLAGS_RELEASE)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-fastsse" GMXC_CXXFLAGS_RELEASE)
    endif()

    # Pathscale
    if (CMAKE_C_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math" 
                         GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_LANG "-std=gnu99" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math" 
                         GMXC_CXXFLAGS_RELEASE)
    endif()

    # xlc
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        GMX_TEST_CFLAG(CFLAGS_OPT "-qarch=auto -qtune=auto" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_LANG "-qlanglvl=extc99" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "XL")
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-qarch=auto -qtune=auto" GMXC_CXXFLAGS)
    endif()

    #msvc
    if (MSVC)
        # disable warnings for: 
        #      inconsistent dll linkage
        GMX_TEST_CFLAG(CFLAGS_WARN "/wd4273" GMXC_CFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4273" GMXC_CXXFLAGS)
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
    endif()

    # Now process nonbonded accelerated kernels settings

    string(TOUPPER ${GMX_CPU_ACCELERATION} ${GMX_CPU_ACCELERATION})
    if(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
        # nothing to do
    elseif(${GMX_CPU_ACCELERATION} STREQUAL "SSE2")

        GMX_TEST_CFLAG(ACCELERATION_SSE2_C_FLAGS "-msse2" ACCELERATION_C_FLAGS)
        if(NOT ACCELERATION_SSE2_C_FLAGS AND GMX_NATIVE_WINDOWS)
            GMX_TEST_CFLAG(ACCELERATION_SSE2_C_FLAGS "/arch:SSE2" ACCELERATION_C_FLAGS)
        endif()

        if (CMAKE_CXX_COMPILER_LOADED)
            GMX_TEST_CXXFLAG(ACCELERATION_SSE2_CXX_FLAGS "-msse2" ACCELERATION_CXX_FLAGS)
            if(NOT ACCELERATION_SSE2_CXX_FLAGS AND GMX_NATIVE_WINDOWS)
                GMX_TEST_CXXFLAG(ACCELERATION_SSE2_CXX_FLAGS "/arch:SSE2" ACCELERATION_CXX_FLAGS)
            endif()
        endif()

        # We don't warn for lacking SSE2 flag support, since that is probably standard today.

        # Only test the include after we have tried to add the correct flag for SSE2 support
        check_include_files(emmintrin.h  HAVE_EMMINTRIN_H ${ACCELERATION_SSE2_C_FLAGS})

        if(NOT HAVE_EMMINTRIN_H)
            message(FATAL_ERROR "Cannot find emmintrin.h, which is required for SSE2 intrinsics support.")
        endif()

        set(GMX_CPU_ACCELERATION_X86_SSE2 1)
        # The user should not be able to set this orthogonally to the acceleration
        set(GMX_X86_SSE2 1)
        if (NOT ACCELERATION_QUIETLY)
            message(STATUS "Enabling SSE2 Gromacs acceleration, and it will help compiler optimization.")
        endif()

    elseif(${GMX_CPU_ACCELERATION} STREQUAL "SSE4.1")

        GMX_TEST_CFLAG(ACCELERATION_SSE4.1_C_FLAGS "-msse4.1" ACCELERATION_C_FLAGS)
        if (NOT ACCELERATION_SSE4.1_C_FLAGS AND GMX_NATIVE_WINDOWS)
            GMX_TEST_CFLAG(ACCELERATION_SSE4.1_C_FLAGS "/arch:SSE4.1" ACCELERATION_C_FLAGS)
        endif()
        if (NOT ACCELERATION_SSE4.1_C_FLAGS)
	    if(MSVC)
                # Not surprising if we end up here! MSVC current does
                # not support the SSE4.1 flag. However, it appears to
                # accept SSE4.1 intrinsics when SSE2 support is
                # enabled, so we try that instead.
                message(WARNING "No C SSE4.1 flag found, but MSVC seems to let GROMACS generate SSE4.1 code if the SSE2 compiler flag is used, so we'll try that.")
                GMX_TEST_CFLAG(ACCELERATION_SSE4.1_C_FLAGS "/arch:SSE2" ACCELERATION_C_FLAGS)
            else()
                message(WARNING "No C SSE4.1 flag found. Consider a newer compiler, or use SSE2 for slightly lower performance.")
            endif()
        endif()

        if (CMAKE_CXX_COMPILER_LOADED)
            GMX_TEST_CXXFLAG(ACCELERATION_SSE4.1_CXX_FLAGS "-msse4.1" ACCELERATION_CXX_FLAGS)
            if (NOT ACCELERATION_SSE4.1_CXX_FLAGS AND GMX_NATIVE_WINDOWS)
                GMX_TEST_CXXFLAG(ACCELERATION_SSE4.1_CXX_FLAGS "/arch:SSE4.1" ACCELERATION_CXX_FLAGS)
            endif()
            if (NOT ACCELERATION_SSE4.1_CXX_FLAGS)
	        if(MSVC)
                    # Not surprising if we end up here! MSVC current does
                    # not support the SSE4.1 flag. However, it appears to
                    # accept SSE4.1 intrinsics when SSE2 support is
                    # enabled, so we try that instead.
                    message(WARNING "No C++ SSE4.1 flag found, but MSVC seems to let GROMACS generate SSE4.1 code if the SSE2 compiler flag is used, so we'll try that.")
                    GMX_TEST_CXXFLAG(ACCELERATION_SSE4.1_CXX_FLAGS "/arch:SSE2" ACCELERATION_CXX_FLAGS)
                else()
                    message(WARNING "No C++ SSE4.1 flag found. Consider a newer compiler, or use SSE2 for slightly lower performance.")
                endif()
            endif()
        endif()

        # This must come after we have added the -msse4.1 flag on some platforms.
        check_include_file(smmintrin.h  HAVE_SMMINTRIN_H ${ACCELERATION_SSE4.1_C_FLAGS})

        if(NOT HAVE_SMMINTRIN_H)
            message(FATAL_ERROR "Cannot find smmintrin.h, which is required for SSE4.1 intrinsics support.")
        endif()

        set(GMX_CPU_ACCELERATION_X86_SSE4_1 1)
        # The user should not be able to set this orthogonally to the acceleration
        set(GMX_X86_SSE4_1 1)
        set(GMX_X86_SSE2   1)
        if (NOT ACCELERATION_QUIETLY)
            message(STATUS "Enabling SSE4.1 Gromacs acceleration, and it will help compiler optimization.")
        endif()

    elseif(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" OR ${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")

        # Set the AVX compiler flag for both these choices!

        GMX_TEST_CFLAG(ACCELERATION_AVX_C_FLAGS "-mavx" ACCELERATION_C_FLAGS)
        if(NOT ACCELERATION_AVX_C_FLAGS AND GMX_NATIVE_WINDOWS)
            GMX_TEST_CFLAG(ACCELERATION_AVX_C_FLAGS "/arch:AVX" ACCELERATION_C_FLAGS)
        endif()
        if(NOT ACCELERATION_AVX_C_FLAGS)
            message(WARNING "No C AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
        endif()

        if (CMAKE_CXX_COMPILER_LOADED)
            GMX_TEST_CXXFLAG(ACCELERATION_AVX_COMMON_CXX_FLAGS "-mavx" ACCELERATION_CXX_FLAGS)
            if(NOT ACCELERATION_AVX_COMMON_CXX_FLAGS AND GMX_NATIVE_WINDOWS)
                GMX_TEST_CXXFLAG(ACCELERATION_AVX_COMMON_CXX_FLAGS "/arch:AVX" ACCELERATION_CXX_FLAGS)
            endif()
            if(NOT ACCELERATION_AVX_COMMON_CXX_FLAGS)
                message(WARNING "No C++ AVX flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
            endif()
        endif()

        # Set the FMA4 flags (MSVC doesn't require any)
        if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA" AND NOT MSVC)
            GMX_TEST_CFLAG(ACCELERATION_AVX_128_FMA_C_FLAGS "-mfma4" ACCELERATION_C_FLAGS)
            if (NOT ACCELERATION_AVX_128_FMA_C_FLAGS)
                message(WARNING "No C FMA4 flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
            endif()

            GMX_TEST_CFLAG(ACCELERATION_AVX_128_XOP_C_FLAGS "-mxop" ACCELERATION_C_FLAGS)
            # No big deal if we do not have xop, so no point yelling warnings about it.

            if (CMAKE_CXX_COMPILER_LOADED)
                GMX_TEST_CXXFLAG(ACCELERATION_AVX_128_FMA_CXX_FLAGS "-mfma4" ACCELERATION_CXX_FLAGS)
                if (NOT ACCELERATION_AVX_128_FMA_CXX_FLAGS)
                    message(WARNING "No C++ FMA4 flag found. Consider a newer compiler, or try SSE4.1 (lower performance).")
                endif()

                GMX_TEST_CXXFLAG(ACCELERATION_AVX_128_XOP_CXX_FLAGS "-mxop" ACCELERATION_CXX_FLAGS)
            endif()
        endif()

        # Only test the header after we have tried to add the flag for AVX support
        check_include_file(immintrin.h  HAVE_IMMINTRIN_H ${ACCELERATION_AVX_C_FLAGS})

        if(NOT HAVE_IMMINTRIN_H)
            message(FATAL_ERROR "Cannot find immintrin.h, which is required for AVX intrinsics support. Consider switching compiler.")
        endif(NOT HAVE_IMMINTRIN_H)

        if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_256")
            try_compile(TEST_AVX ${CMAKE_BINARY_DIR}
                "${CMAKE_SOURCE_DIR}/cmake/TestAVX.c"
                COMPILE_DEFINITIONS "${ACCELERATION_AVX_C_FLAGS}")
            if(NOT TEST_AVX)
                message(FATAL_ERROR "Cannot compile AVX intrinsics. Consider switching compiler.")
            endif()
        endif()

        # GCC requires x86intrin.h for FMA support. MSVC 2010 requires intrin.h for FMA support.
        check_include_file(x86intrin.h HAVE_X86INTRIN_H ${ACCELERATION_AVX_C_FLAGS})
        check_include_file(intrin.h HAVE_INTRIN_H ${ACCELERATION_AVX_C_FLAGS})

        # The user should not be able to set this orthogonally to the acceleration
        set(GMX_X86_SSE4_1 1)
        set(GMX_X86_SSE2   1)

        # But just enable one of the choices internally...
        if(${GMX_CPU_ACCELERATION} STREQUAL "AVX_128_FMA")
            set(GMX_CPU_ACCELERATION_X86_AVX_128_FMA 1)
            set(GMX_X86_AVX_128_FMA 1)
            if (NOT ACCELERATION_QUIETLY)
                message(STATUS "Enabling 128-bit AVX Gromacs acceleration (with fused-multiply add), and it will help compiler optimization.")
            endif()
        else()
            # If we are not doing AVX_128, it must be AVX_256...
            set(GMX_CPU_ACCELERATION_X86_AVX_256 1)
            set(GMX_X86_AVX_256 1)
            if (NOT ACCELERATION_QUIETLY)
                message(STATUS "Enabling 256-bit AVX Gromacs acceleration, and it will help compiler optimization.")
            endif()
        endif()

    elseif(${GMX_CPU_ACCELERATION} STREQUAL "FORTRAN")

        #    Fortran is temporarily disabled while we push in nbNxN kernels.
        #    We need to fake it a bit here to avoid jenkins build errors!
        #    add_definitions(-DGMX_FORTRAN)

    elseif(${GMX_CPU_ACCELERATION} STREQUAL "BLUEGENE")
        # GMX_CPU_ACCELERATION=BlueGene should be set in the Toolchain-BlueGene?-???.cmake file
        if (NOT ACCELERATION_QUIETLY)
            message(STATUS "Configuring for BlueGene")
        endif()
        set(GMX_BLUEGENE 1)
        if (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
            set(SHARED_LIBS_DEFAULT OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
            set(BUILD_SHARED_LIBS OFF CACHE BOOL "Shared libraries not compatible with BlueGene/L, disabled!" FORCE)
        endif (${CMAKE_SYSTEM_NAME} STREQUAL "BlueGeneL")
        set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on BlueGene" FORCE)
        set(GMX_POWERPC_INVSQRT ON CACHE BOOL "Use hardware reciprocal square root on BlueGene" FORCE)
        set(GMX_X11 OFF CACHE BOOL "X11 not compatible with BlueGene, disabled!" FORCE)
        set(GMX_THREAD_MPI OFF CACHE BOOL "Thread-MPI not compatible with BlueGene, disabled!" FORCE)
        set(GMX_MPI ON CACHE BOOL "Use MPI on BlueGene" FORCE)
        # Access to /etc/passwd is not available on the back end of BlueGene,
        # despite being detected by CMake. This can cause linker warnings
        # about harmless things in src/gmxlib/string2.h.
        set(HAVE_PWD_H OFF)
        # The automatic testing for endianness does not work for the BlueGene cross-compiler
        set(GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER 1 CACHE INTERNAL "BlueGene has big endian FP byte order (by default)" FORCE)
        set(GMX_IEEE754_BIG_ENDIAN_WORD_ORDER 1 CACHE INTERNAL "BlueGene has big endian FP word order (by default)" FORCE)
    elseif(${GMX_CPU_ACCELERATION} STREQUAL "POWER6")
        set(GMX_POWER6 1)
        set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Do not use software reciprocal square root on Power6" FORCE)
        set(GMX_POWERPC_INVSQRT ON CACHE BOOL "Use hardware reciprocal square root on Power6" FORCE)
    else(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
        MESSAGE(FATAL_ERROR "Unrecognized option for accelerated kernels: ${GMX_CPU_ACCELERATION}. Pick one of None, SSE2, SSE4.1, AVX_128_FMA, AVX_256, Fortran, BlueGene, Power6")
    endif(${GMX_CPU_ACCELERATION} STREQUAL "NONE")
    set(ACCELERATION_QUIETLY TRUE CACHE INTERNAL "")

    ########################################################################
    # Fix stupid flags on Windows
    ########################################################################
    IF( WIN32 AND NOT CYGWIN)
        if (NOT BUILD_SHARED_LIBS)
            option(GMX_PREFER_STATIC_LIBS "When finding libraries prefer static system libraries (MT instead of MD)!" ON)
            mark_as_advanced(GMX_PREFER_STATIC_LIBS)
            SET(SHARED_LIBS_DEFAULT OFF)
        else()
            add_definitions(-DUSE_VISIBILITY -DTMPI_USE_VISIBILITY)
        endif()

        IF (GMX_PREFER_STATIC_LIBS)
            #Only setting Debug and Release flags. Others configurations current not used.
            STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
            SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
            STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
            SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)
            if(CMAKE_CXX_COMPILER_LOADED)
                STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
                SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} CACHE STRING "" FORCE)
                STRING(REPLACE /MD /MT CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
                SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} CACHE STRING "" FORCE)
            endif()
        ENDIF()

        #Workaround for cmake bug 13174. Replace deprecated options.
        IF( CMAKE_C_COMPILER_ID MATCHES "Intel" )
            STRING(REPLACE /GZ /RTC1 CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
            SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)
        ENDIF()
        IF( CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_LOADED)
            STRING(REPLACE /GZ /RTC1 CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
            STRING(REPLACE /GX /EHsc CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
            SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} CACHE STRING "" FORCE)

            STRING(REPLACE /GX /EHsc CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
            SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} CACHE STRING "" FORCE)
        ENDIF()
    ENDIF()

ENDMACRO(gmx_c_flags)

