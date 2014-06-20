#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

# Test C flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CFLAGSVAR.
MACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_C_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF()
    IF (${VARIABLE})
        SET (${CFLAGSVAR} "${FLAGS} ${${CFLAGSVAR}}")
    ENDIF ()
ENDMACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)

# Test C++ flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CXXFLAGSVAR.
MACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_CXX_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF()
    IF (${VARIABLE})
        SET (${CXXFLAGSVAR} "${FLAGS} ${${CXXFLAGSVAR}}")
    ENDIF ()
ENDMACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)

# Set the real CMake variables for compiler flags. This should be a function
# so we can have proper local variables while avoiding duplicating code.
function(gmx_set_cmake_compiler_flags)
    foreach(language C CXX)
        # Copy the flags for the release build type to the build types
        # that are modified forms of it. Ideally, the list of build
        # types that are modifications of the Release build type would
        # be set up elsewhere and passed to this function, but it is
        # inconvenient in CMake to pass more than one list, and such a
        # list is only used here.
        foreach(build_type RELWITHDEBUGINFO RELWITHASSERT MINSIZEREL)
            set(GMXC_${language}FLAGS_${build_type} "${GMXC_${language}FLAGS_RELEASE}")
        endforeach()
        # Copy the flags that are only used by the real Release build
        # type. Currently unused, but we plan to use -Wno-array-bounds
        # in Release to work around gcc-4.8 being a little too vocal
        # about some perfectly good code, while using RelWithAssert
        # (ie. without that suppression) in Jenkins.
        set(GMXC_${language}FLAGS_RELEASE "${GMXC_${language}FLAGS_RELEASE} ${GMXC_${language}FLAGS_RELEASE_ONLY}")

        # Modify the real CMake variables for compiler flags for all
        # builds and language types, and also those common to all
        # build types.
        foreach(build_type "" ${build_types_with_explicit_flags})
            if("${build_type}" STREQUAL "")
                set(punctuation "") # for general compiler flags (e.g.) CMAKE_CXX_FLAGS
            else()
                set(punctuation "_") # for build-type-specific compiler flags (e.g.) CMAKE_CXX_FLAGS_RELEASE
            endif()

            # Append to the variables for the given build type for
            # each language, in the parent scope.
            set(CMAKE_${language}_FLAGS${punctuation}${build_type}
                "${GMXC_${language}FLAGS${punctuation}${build_type}} ${CMAKE_${language}_FLAGS${punctuation}${build_type}}"
                PARENT_SCOPE)
        endforeach()
    endforeach()
endfunction()

# This is the actual exported function to be called 
MACRO(gmx_c_flags)

    include(CheckCCompilerFlag)
    include(CheckCXXCompilerFlag)

    # gcc
    if(CMAKE_COMPILER_IS_GNUCC)
        #flags are added in reverse order and -Wno* need to appear after -Wall
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value -Wunused-parameter" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare -Wpointer-arith" GMXC_CFLAGS)
        # Since 4.8 on by default. For previous version disabling is a no-op. Only disabling for Release because with assert
        # the warnings are OK.
        GMX_TEST_CFLAG(CFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CFLAGS_RELEASE_ONLY)
        # Since gcc 4.8 strict - false postives with old gmx_fatal. TODO: Remove in master
        GMX_TEST_CFLAG(CFLAGS_WARN_UNINIT "-Wno-maybe-uninitialized" GMXC_CFLAGS)
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
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused-function" GMXC_CXXFLAGS)
        # Problematic with CUDA
        # GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EFFCXX "-Wnon-virtual-dtor" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wpointer-arith" GMXC_CXXFLAGS)
        GMX_TEST_CFLAG(CXXFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CXXFLAGS_RELEASE_ONLY)
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
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "-wd161" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_WARN "-w3 -wd111 -wd177 -wd181 -wd193 -wd271 -wd304 -wd383 -wd424 -wd444 -wd522 -wd593 -wd869 -wd981 -wd1418 -wd1419 -wd1572 -wd1599 -wd2259 -wd2415 -wd2547 -wd2557 -wd3280 -wd3346" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_STDGNU "-std=gnu99" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_OPT "-ip -funroll-all-loops -alias-const -ansi-alias" GMXC_CFLAGS_RELEASE)
        else()
            if(NOT GMX_OPENMP)
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "/wd161" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_WARN "/W3 /wd111 /wd177 /wd181 /wd193 /wd271 /wd304 /wd383 /wd424 /wd444 /wd522 /wd593 /wd869 /wd981 /wd1418 /wd1419 /wd1572 /wd1599 /wd2259 /wd2415 /wd2547 /wd2557 /wd3280 /wd3346" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_OPT "/Qip" GMXC_CFLAGS_RELEASE)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-wd161" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-w3 -wd111 -wd177 -wd181 -wd193 -wd271 -wd304 -wd383 -wd424 -wd444 -wd522 -wd593 -wd869 -wd981 -wd1418 -wd1419 -wd1572 -wd1599 -wd2259 -wd2415 -wd2547 -wd2557 -wd3280 -wd3346 -wd1782" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-ip -funroll-all-loops -alias-const -ansi-alias" GMXC_CXXFLAGS_RELEASE)
        else()
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "/wd161" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/W3 /wd111 /wd177 /wd181 /wd193 /wd271 /wd304 /wd383 /wd424 /wd444 /wd522 /wd593 /wd869 /wd981 /wd1418 /wd1419 /wd1572 /wd1599 /wd2259 /wd2415 /wd2547 /wd2557 /wd3280 /wd3346 /wd1782" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "/Qip" GMXC_CXXFLAGS_RELEASE)
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

    # msvc
    if (MSVC)
        # disable warnings for: 
        #      forcing value to bool
        #      "this" in initializer list
        #      deprecated (posix, secure) functions
        #      truncation (double -> float)
        #      conversion from 'double' to 'real', possible loss of data
        #      unreferenced local variable (only C)
        #      conversion from 'size_t' to 'int', possible loss of data
        #      conversion from 'const char*' to 'void*', different 'const' qualifiers (only C)
        GMX_TEST_CFLAG(CFLAGS_WARN "/wd4800 /wd4355 /wd4996 /wd4305 /wd4244 /wd4101 /wd4267 /wd4090" GMXC_CFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4800 /wd4355 /wd4996 /wd4305 /wd4244 /wd4267" GMXC_CXXFLAGS)
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value -Wunused-parameter" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wpointer-arith" GMXC_CFLAGS_EXTRA)
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused-function" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wpointer-arith" GMXC_CXXFLAGS)
    endif()

    # now actually set the flags:
    if (NOT GMX_SKIP_DEFAULT_CFLAGS)
        gmx_set_cmake_compiler_flags()
    endif()
ENDMACRO(gmx_c_flags)

