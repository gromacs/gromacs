#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
        foreach(build_type RELWITHDEBINFO RELWITHASSERT MINSIZEREL PROFILE)
            set(GMXC_${language}FLAGS_${build_type} "${GMXC_${language}FLAGS_RELEASE}")
        endforeach()
        # Copy the flags that are only used by the real Release build
        # type. Used for, e.g., -Wno-array-bounds in Release to work around
        # gcc-4.8 being a little too vocal about some perfectly good code,
        # while using RelWithAssert (ie. without that suppression) in Jenkins.
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
            # each language, in the parent scope. We add our new variables at the end, so
            # compiler-specific choices are more likely to override default CMake choices.
            # This is for instance useful for RelWithDebInfo builds, where we want to use the full
            # set of our optimization flags detected in this file, rather than having -O2 override them.
            set(CMAKE_${language}_FLAGS${punctuation}${build_type}
                "${CMAKE_${language}_FLAGS${punctuation}${build_type}} ${GMXC_${language}FLAGS${punctuation}${build_type}}"
                PARENT_SCOPE)
        endforeach()
    endforeach()
endfunction()

# This is the actual exported function to be called
macro (gmx_c_flags)

    include(CheckCCompilerFlag)
    include(CheckCXXCompilerFlag)

    # gcc
    if(CMAKE_COMPILER_IS_GNUCC)
        #flags are added in reverse order and -Wno* need to appear after -Wall
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value -Wunused-parameter" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare -Wpointer-arith" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_UNDEF "-Wundef" GMXC_CFLAGS)
            # Since 4.8 on by default. For previous version disabling is a no-op. Only disabling for Release because with assert
            # the warnings are OK.
            GMX_TEST_CFLAG(CFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CFLAGS_RELEASE_ONLY)
            if(CYGWIN)
                GMX_TEST_CFLAG(CFLAGS_WARN_SUBSCRIPT "-Wno-char-subscripts" GMXC_CFLAGS)
            endif()
        endif()
        # new in gcc 4.5
        GMX_TEST_CFLAG(CFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_COPT "-funroll-all-loops"
                       GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_NOINLINE "-fno-inline" GMXC_CFLAGS_DEBUG)
    endif()
    # g++
    if(CMAKE_COMPILER_IS_GNUCXX)
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            # Problematic with CUDA
            # GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EFFCXX "-Wnon-virtual-dtor" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wpointer-arith -Wmissing-declarations" GMXC_CXXFLAGS)
            # CUDA versions prior to 7.5 come with a header (math_functions.h) which uses the _MSC_VER macro
            # unconditionally, so we don't use -Wundef for earlier CUDA versions.
            if(NOT(GMX_GPU AND CUDA_VERSION VERSION_LESS "7.5"))
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN_UNDEF "-Wundef" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CFLAG(CXXFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CXXFLAGS_RELEASE_ONLY)
        endif()
        # new in gcc 4.5
        GMX_TEST_CXXFLAG(CXXFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_COPT "-funroll-all-loops"
                         GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_NOINLINE "-fno-inline" GMXC_CXXFLAGS_DEBUG)
    endif()

    # icc
    if (CMAKE_C_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32)
            if(NOT GMX_OPENMP)
# 3180: unrecognized OpenMP #pragma
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "-wd3180" GMXC_CFLAGS)
            endif()
            if (GMX_COMPILER_WARNINGS)
# -w3 enables a lot of useful diagnostics but we don't care about all. -wd disables some selectively.
# 177: function/variable ".." was declared but never referenced
# 411: class defines no constructor to initialize the following (incorrect for struct, initializer list works)
# 593: variable ".." was set but never used
# 981: operands are evaluated in unspecified order
#1418: external function definition with no prior declaration
#1419: external declaration in primary source file
#1572: floating-point equality and inequality comparisons are unreliable
#1599: declaration hides variable ".."
#2259: non-pointer conversion from ".." to ".." may lose significant bits
#2415: variable ".." of static storage duration was declared but never referenced
#2547: ".." was specified as both a system and non-system include directory
#2557: comparison between signed and unsigned operands
#3280: declaration hides member ".."
#11074: Inlining inhibited by limit max-size(/max-total-size)
#11076: To get full report use -opt-report=3 -opt-report-phase ipo (shown for previous remark)
                GMX_TEST_CFLAG(CFLAGS_WARN "-w3 -wd177 -wd411 -wd593 -wd981 -wd1418 -wd1419 -wd1572 -wd1599 -wd2259 -wd2415 -wd2547 -wd2557 -wd3280 -wd11074 -wd11076" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_STDGNU "-std=gnu99" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_OPT "-ip -funroll-all-loops -alias-const -ansi-alias -no-prec-div -fimf-domain-exclusion=14" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_DEBUG "-O0" GMXC_CFLAGS_DEBUG) #icc defaults to -O2 even with -g
            GMX_TEST_CFLAG(CFLAGS_FP_RELASSERT "-fp-model except -fp-model precise" GMXC_CFLAGS_RELWITHASSERT)
        else()
            if(NOT GMX_OPENMP)
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "/wd3180" GMXC_CFLAGS)
            endif()
            if (GMX_COMPILER_WARNINGS)
GMX_TEST_CFLAG(CFLAGS_WARN "/W3 /wd177 /wd411 /wd593 /wd981 /wd1418 /wd1419 /wd1572 /wd1599 /wd2259 /wd2415 /wd2547 /wd2557 /wd3280" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_OPT "/Qip" GMXC_CFLAGS_RELEASE)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-wd3180" GMXC_CXXFLAGS)
            endif()
            if (GMX_COMPILER_WARNINGS)
                if (GMX_GPU)
# Suppress warnings from CUDA headers
# 7:   unrecognized token
# 82:  storage class is not first
# The below are also required for math_functions.h / math_functions.hpp at least until CUDA 8.0-RC
# 193: zero used for undefined preprocessing identifer
# 3346:dynamic exception specifiers are deprecated
                    GMX_TEST_CXXFLAG(CXXFLAGS_WARN_OLD_GPU "-wd7 -wd82 -wd193 -wd3346" GMXC_CXXFLAGS)
                endif()
#All but the following warnings are identical for the C-compiler (see above)
# 383: value copied to temporary, reference to temporary used
# 444: destructor for base class ".." is not virtual
#2282: unrecognized GCC pragma
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-w3 -wd177 -wd383 -wd411 -wd444 -wd981 -wd1418 -wd1572 -wd1599 -wd2259 -wd2547 -wd3280 -wd11074 -wd11076 -wd2282" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-ip -funroll-all-loops -alias-const -ansi-alias -no-prec-div -fimf-domain-exclusion=14" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_DEBUG "-O0" GMXC_CXXFLAGS_DEBUG)
            GMX_TEST_CXXFLAG(CXXFLAGS_FP_RELASSERT "-fp-model except -fp-model precise" GMXC_CXXFLAGS_RELWITHASSERT)
        else()
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "/wd3180" GMXC_CFLAGS)
            endif()
            if (GMX_COMPILER_WARNINGS)
#809: exception specification for virtual function X is incompatible with that of overridden function
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/W3 /wd177 /wd383 /wd411 /wd444 /wd809 /wd981 /wd1418 /wd1572 /wd1599 /wd2259 /wd2547 /wd3280 /wd11074 /wd11076 /wd2282" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "/Qip" GMXC_CXXFLAGS_RELEASE)
        endif()
    endif()

    # PGI
    # Inter-procedural analysis causes pgcc/pgc++ to crash when linking the library with PGI release 15.7.
    if (CMAKE_C_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CFLAG(CFLAGS_OPT "-Mnoipa" GMXC_CFLAGS_RELEASE)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
        # Using ipa exposes internal PGI-15.7 compiler bugs at compile time
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-Mnoipa" GMXC_CXXFLAGS_RELEASE)
        # PGI identifies itself as GCC, but does not understand the GCC
        # pragmas that occur in parser.cpp. Since that file is generated
        # we cannot add a define there, but supress the warning instead.
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "--diag_suppress=1675" GMXC_CXXFLAGS)
    endif()

    # Pathscale
    if (CMAKE_C_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math"
                         GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_LANG "-std=gnu99" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math"
                         GMXC_CXXFLAGS_RELEASE)
    endif()

    # xlc
    # The suppressions below stop
    # 1500-036: (I) about -O3 causing non-strict IEEE compliance that changes the semantics of the program (duh)
    # 1500-010: (W) about correct PBC-related use of maximum array indices of DIM-sized C arrays
    # 1500-030: (I) Additional optimization may be attained by recompiling and specifying MAXMEM option with a value greater than 8192.
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        GMX_TEST_CFLAG(CFLAGS_ARCH "-qarch=auto -qtune=auto" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_OPT  "-O3" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_LANG "-qlanglvl=extc99" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_LANG "-qsuppress=1500-036 -qsuppress=1500-010 -qsuppress=1500-030" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "XL")
        GMX_TEST_CXXFLAG(CXXFLAGS_ARCH "-qarch=auto -qtune=auto" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-O3" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_LANG "-qsuppress=1500-036 -qsuppress=1500-010 -qsuppress=1500-030" GMXC_CXXFLAGS)
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
        if(NOT CMAKE_CONFIGURATION_TYPES)
            GMX_TEST_CFLAG(CFLAGS_WARN "/wd4800 /wd4355 /wd4996 /wd4305 /wd4244 /wd4101 /wd4267 /wd4090" GMXC_CFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4800 /wd4355 /wd4996 /wd4305 /wd4244 /wd4267" GMXC_CXXFLAGS)
        else() #Projects only use the C++ flags
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4800 /wd4355 /wd4996 /wd4305 /wd4244 /wd4101 /wd4267 /wd4090" GMXC_CXXFLAGS)
        endif()
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value -Wunused-parameter" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wpointer-arith" GMXC_CFLAGS_EXTRA)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wpointer-arith -Wmissing-prototypes" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_DEPRECATED "-Wdeprecated" GMXC_CXXFLAGS)
        endif()
    endif()

    # Fujitsu compilers on PrimeHPC/Sparc64
    if(${CMAKE_C_COMPILER_ID} MATCHES Fujitsu OR
       (${CMAKE_C_COMPILER_ID} MATCHES unknown AND ${CMAKE_C_COMPILER} MATCHES ^fcc))
        GMX_TEST_CFLAG(CFLAG_GNUCOMPAT "-Xg -w" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAG_OPT "-Kfast,reduction,swp,simd=2,uxsimd,fsimple -x100" GMXC_CFLAGS)
    endif()

    if(${CMAKE_CXX_COMPILER_ID} MATCHES Fujitsu OR
       (${CMAKE_CXX_COMPILER_ID} MATCHES unknown AND ${CMAKE_CXX_COMPILER} MATCHES ^FCC))
        GMX_TEST_CXXFLAG(CXXFLAG_GNUCOMPAT "-Xg -w" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAG_OPT "-Kfast,reduction,swp,simd=2,uxsimd,fsimple -x100" GMXC_CXXFLAGS)
    endif()

    # now actually set the flags:
    if (NOT GMX_SKIP_DEFAULT_CFLAGS)
        gmx_set_cmake_compiler_flags()
    endif()

endmacro()
