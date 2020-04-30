#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013 by the GROMACS development team.
# Copyright (c) 2014,2015,2016,2017,2018 by the GROMACS development team.
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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
        list(APPEND ${CFLAGSVAR} "${FLAGS}")
    ENDIF ()
ENDMACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)

# Test C++ flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CXXFLAGSVAR.
MACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_CXX_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF()
    IF (${VARIABLE})
        list(APPEND ${CXXFLAGSVAR} "${FLAGS}")
    ENDIF ()
ENDMACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)

# Prepare some local variables so CUDA and non-CUDA code in targets
# works the same way.
function(gmx_target_compile_options_inner)
    set (CFLAGS "${SIMD_C_FLAGS};${MPI_COMPILE_FLAGS};${EXTRA_C_FLAGS};${GMXC_CFLAGS}" PARENT_SCOPE)
    set (CXXFLAGS "${SIMD_CXX_FLAGS};${MPI_COMPILE_FLAGS};${EXTRA_CXX_FLAGS};${GMXC_CXXFLAGS}" PARENT_SCOPE)
endfunction()

# Implementation function to add compiler flags expected for all
# GROMACS build configurations, and those expected for the current
# CMake build type (e.g. Release) to TARGET. Other GROMACS CMake code
# is expected to use either gmx_target_compile_options(name_of_target)
# or gmx_cuda_target_compile_options(name_of_variable) because CUDA
# compilation has special requirements.
#
# Most targets (ie. libraries, executables) need compiler flags that
# are characteristic of the build configuration. This function
# provides a central point for adding such flags. This approach is
# more flexible than e.g. setting CMAKE_CXX_FLAGS globally, because
# that setting will apply globally, which means it applies also to
# "external" code that the build of GROMACS might also build.
function(gmx_target_compile_options TARGET)
    if (GMX_SKIP_DEFAULT_CFLAGS)
        return()
    endif()

    # Prepare the generic compiler options
    gmx_target_compile_options_inner()
    target_compile_options(${TARGET}
        PRIVATE
        $<$<COMPILE_LANGUAGE:C>:${CFLAGS}>
        $<$<COMPILE_LANGUAGE:CXX>:${CXXFLAGS}>
        )
    # Add compiler options for the build types
    foreach(build_type ${build_types_with_explicit_flags})
        target_compile_options(${TARGET}
            BEFORE PRIVATE
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:${build_type}>>:${GMXC_CFLAGS_${build_type}}>
            $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:${build_type}>>:${GMXC_CXXFLAGS_${build_type}}>
            )
    endforeach()
    # Add the release-configuration compiler options to build
    # configurations that derive from it.
    foreach(build_type RELWITHDEBINFO RELWITHASSERT MINSIZEREL PROFILE)
        target_compile_options(${TARGET}
            BEFORE PRIVATE
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:${build_type}>>:${GMXC_CFLAGS_RELEASE}>
            $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:${build_type}>>:${GMXC_CXXFLAGS_RELEASE}>
            )
    endforeach()
    # Add those flags that we only want in the proper release build
    # configuration.
    target_compile_options(${TARGET}
        BEFORE PRIVATE
        $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:RELEASE>>:${GMXC_CFLAGS_RELEASE_ONLY}>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:${GMXC_CXXFLAGS_RELEASE_ONLY}>
        )
endfunction()

# The approach taken by FindCUDA.cmake is to require that the compiler
# flags are known and present in a variable before creating the target
# for the library. (Those get embedded in files that are generated at
# the point of calling cuda_add_library, which does not create a
# target that one could later call target_compile_options upon.) So,
# this function instead returns appropriate content in
# ${VARIABLE_NAME}, along with other such variables that are
# specialized for the various build_types. Hopefully this will improve
# when we use native CUDA language support in our CMake.
function(gmx_cuda_target_compile_options VARIABLE_NAME)
    if (GMX_SKIP_DEFAULT_CFLAGS)
        set (CXXFLAGS "")
    else()
        # Prepare the generic compiler options
        gmx_target_compile_options_inner()
        # CUDA headers issue lots of warnings when compiled with
        # -Wundef because they use old-style #ifdef a lot. We'd prefer
        # to have FindCUDA.cmake treat CUDA internal headers with
        # -isystem so that these warnings are naturally suppressed,
        # but there's no way to do that without bundling a modified
        # form of FindCUDA.cmake. That creates its own problems,
        # because people either don't know we do that, or don't
        # remember that we don't do that in user tarballs.
        #
        # We have make check-source ensuring that we have included
        # config.h any time we use such symbols in commits in a merge
        # request. Local development could run that too. So, we can
        # tolerate any remaining risk from accidentally using
        # e.g. #ifdef GMX_MPI rather than #if GMX_MPI in CUDA source
        # files.
        #
        # So we disable -Wundef by the simple hack of appending
        # -Wno-undef after it. That's more maintainable than having
        # logic to avoid adding -Wundef to GMXC_CXXFLAGS, given the
        # current approach to adding them. Hopefully this will improve
        # if/when we have more CMake object libraries, and/or native
        # CUDA compilation.
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NOUNDEF "-Wno-undef" CXXFLAGS)
    endif()

    # Only C++ compilation is supported with CUDA code in GROMACS
    set(${VARIABLE_NAME} ${CXXFLAGS} PARENT_SCOPE)

    # Now organize the flags for different build
    # configurations. First, the debug configuration.
    set(${VARIABLE_NAME}_DEBUG "${GMXC_CXXFLAGS_DEBUG}" PARENT_SCOPE)
    # Add those flags that we only want in the proper release build
    # configuration.
    set(${VARIABLE_NAME}_RELEASE "${GMXC_CXXFLAGS_RELEASE};${GMXC_CXXFLAGS_RELEASE_ONLY}" PARENT_SCOPE)
    # Add the release flags to build configurations that derive from it
    foreach(build_type RELWITHDEBINFO RELWITHASSERT MINSIZEREL PROFILE)
        set(${VARIABLE_NAME}_${build_type} "${GMXC_CXXFLAGS_RELEASE};${GMXC_CXXFLAGS_${build_type}}" PARENT_SCOPE)
    endforeach()
endfunction()

# Add the WARNING_FLAG to the compile options for TARGET if the
# compiler supports that flag. VARNAME is used to cache the result of
# the check whether the compiler supports that flag, as multiple
# targets may use the same suppression.
#
# This is generally used to suppress warnings judged not worth fixing
# in code external to, or generated by, GROMACS, or code that is
# more efficient to work around and later replace, rather than fix.
function(gmx_target_warning_suppression TARGET WARNING_FLAG VARNAME)
    check_cxx_compiler_flag(${WARNING_FLAG} ${VARNAME})
    if(${VARNAME})
        target_compile_options(${TARGET} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${WARNING_FLAG}>)
    endif()
endfunction()

# Add the WARNING_FLAG to the compile flags for SOURCE_FILE if the
# compiler supports that flag. VARNAME is used to cache the result of
# the check whether the compiler supports that flag, as multiple
# targets may use the same suppression.
#
# This is generally used to suppress warnings judged not worth fixing
# in code external to, or generated by, GROMACS, or code that is
# more efficient to work around and later replace, rather than fix.
function(gmx_source_file_warning_suppression SOURCE_FILE WARNING_FLAG VARNAME)
    check_cxx_compiler_flag(${WARNING_FLAG} ${VARNAME})
    if(${VARNAME})
        set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS ${WARNING_FLAG})
    endif()
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
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall;-Wno-unused;-Wunused-value;-Wunused-parameter" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wextra;-Wno-missing-field-initializers;-Wno-sign-compare;-Wpointer-arith" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_UNDEF "-Wundef" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CFLAGS_RELEASE_ONLY)
            if(CYGWIN)
                GMX_TEST_CFLAG(CFLAGS_WARN_SUBSCRIPT "-Wno-char-subscripts" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_STRINGOP_TRUNCATION "-Werror=stringop-truncation" GMXC_CFLAGS)
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
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra;-Wno-missing-field-initializers;-Wpointer-arith;-Wmissing-declarations" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_UNDEF "-Wundef" GMXC_CXXFLAGS)
            GMX_TEST_CFLAG(CXXFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CXXFLAGS_RELEASE_ONLY)
            GMX_TEST_CXXFLAG(CXXFLAGS_STRINGOP_TRUNCATION "-Wstringop-truncation" GMXC_CXXFLAGS)
            if (CXXFLAGS_STRINGOP_TRUNCATION)
                set(CXXFLAGS_NO_STRINGOP_TRUNCATION "-Wno-stringop-truncation")
            endif()
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
# 280: selector expression is constant
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
                GMX_TEST_CFLAG(CFLAGS_WARN "-w3;-wd177;-wd280;-wd411;-wd593;-wd981;-wd1418;-wd1419;-wd1572;-wd1599;-wd2259;-wd2415;-wd2547;-wd2557;-wd3280;-wd11074;-wd11076" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_STDGNU "-std=gnu99" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_OPT "-ip;-funroll-all-loops;-alias-const;-ansi-alias;-no-prec-div;-fimf-domain-exclusion=14;-qoverride-limits" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_DEBUG "-O0" GMXC_CFLAGS_DEBUG) #icc defaults to -O2 even with -g
            # The "except" fp-model requires something other than the
            # default "fast" model, so we test and use it with
            # "precise".
            GMX_TEST_CFLAG(CFLAGS_FP_MODEL_RELASSERT "-fp-model=except;-fp-model=precise" GMXC_CFLAGS_RELWITHASSERT)
        else()
            if(NOT GMX_OPENMP)
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "/wd3180" GMXC_CFLAGS)
            endif()
            if (GMX_COMPILER_WARNINGS)
#only on Windows
#161: unrecognized pragma
#1786 function was declared deprecated (is issued for stdlib function such as strncpy which have a _s version)
GMX_TEST_CFLAG(CFLAGS_WARN "/W3;/wd161;/wd177;/wd411;/wd593;/wd981;/wd1418;/wd1419;/wd1572;/wd1599;/wd1786;/wd2259;/wd2415;/wd2547;/wd2557;/wd3280" GMXC_CFLAGS)
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
#All but the following warnings are identical for the C-compiler (see above)
# 304: access control not specified
# 383: value copied to temporary, reference to temporary used
# 444: destructor for base class ".." is not virtual
# 869: was never referenced (false positives)
#2282: unrecognized GCC pragma
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-w3;-wd177;-wd280;-wd304;-wd383;-wd411;-wd444;-wd869;-wd981;-wd1418;-wd1572;-wd1599;-wd2259;-wd2547;-wd3280;-wd11074;-wd11076;-wd2282" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-ip;-funroll-all-loops;-alias-const;-ansi-alias;-no-prec-div;-fimf-domain-exclusion=14;-qoverride-limits" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_DEBUG "-O0" GMXC_CXXFLAGS_DEBUG)
            # The "except" fp-model requires something other than the
            # default "fast" model, so we test and use it with
            # "precise".
            GMX_TEST_CXXFLAG(CXXFLAGS_FP_MODEL_RELASSERT "-fp-model=except;-fp-model=precise" GMXC_CXXFLAGS_RELWITHASSERT)
        else()
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "/wd3180" GMXC_CFLAGS)
            endif()
            if (GMX_COMPILER_WARNINGS)
#161: unrecognized pragma
#809: exception specification for virtual function X is incompatible with that of overridden function
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/W3;/wd161;/wd177;/wd280;/wd304;/wd383;/wd411;/wd444;/wd809;/wd869;/wd981;/wd1418;/wd1572;/wd1599;/wd1786;/wd2259;/wd2547;/wd3280;/wd11074;/wd11076;/wd2282" GMXC_CXXFLAGS)
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
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall;-Wno-unused;-Wunused-value" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_OPT "-OPT:Ofast;-fno-math-errno;-ffast-math"
                         GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_LANG "-std=gnu99" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall;-Wno-unused;-Wunused-value" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-OPT:Ofast;-fno-math-errno;-ffast-math"
                         GMXC_CXXFLAGS_RELEASE)
    endif()

    # xlc
    # The suppressions below stop
    # 1500-036: (I) about -O3 causing non-strict IEEE compliance that changes the semantics of the program (duh)
    # 1500-010: (W) about correct PBC-related use of maximum array indices of DIM-sized C arrays
    # 1500-030: (I) Additional optimization may be attained by recompiling and specifying MAXMEM option with a value greater than 8192.
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        GMX_TEST_CFLAG(CFLAGS_ARCH "-qarch=auto;-qtune=auto" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_OPT  "-O3" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_LANG "-qlanglvl=extc99" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_LANG "-qsuppress=1500-036;-qsuppress=1500-010;-qsuppress=1500-030" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "XL")
        GMX_TEST_CXXFLAG(CXXFLAGS_ARCH "-qarch=auto;-qtune=auto" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-O3" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_LANG "-qsuppress=1500-036;-qsuppress=1500-010;-qsuppress=1500-030" GMXC_CXXFLAGS)
    endif()

    # msvc
    if (MSVC)
        # disable warnings for: 
        #      forcing value to bool
        #      "this" in initializer list
        #      deprecated (posix, secure) functions
        #      C4305: truncation (double -> float)
        #      C4244: conversion from '.*' to '.*', possible loss of data
        #      unreferenced local variable (only C)
        #      C4267: conversion from 'size_t' to 'int', possible loss of data
        #      conversion from 'const char*' to 'void*', different 'const' qualifiers (only C)
        #      unknown pragma (4068)
        #      remark #280: selector expression is constant
        if(NOT CMAKE_CONFIGURATION_TYPES)
            GMX_TEST_CFLAG(CFLAGS_WARN "/wd4800;/wd4355;/wd4996;/wd4305;/wd4244;/wd4101;/wd4267;/wd4090;/wd4068;" GMXC_CFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4800;/wd4355;/wd4996;/wd4305;/wd4244;/wd4267;/wd4068;" GMXC_CXXFLAGS)
        else() # MSVC projects only use the C++ flags
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4800;/wd4355;/wd4996;/wd4305;/wd4244;/wd4101;/wd4267;/wd4090;/wd4068;" GMXC_CXXFLAGS)
        endif()
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall;-Wno-unused;-Wunused-value;-Wunused-parameter" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wpointer-arith" GMXC_CFLAGS_EXTRA)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if (GMX_COMPILER_WARNINGS)
            # If used, -Wall should precede other options that silence warnings it enables
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.0") #LLVM BUG #21629
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_BRACES "-Wno-missing-braces" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra;-Wno-missing-field-initializers;-Wpointer-arith;-Wmissing-prototypes" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_DEPRECATED "-Wdeprecated" GMXC_CXXFLAGS)
            # Functions placed in headers for inlining are not always
            # used in every translation unit that includes the files,
            # so we must disable the warning that there are such
            # functions that are unused.
            GMX_TEST_CXXFLAG(CXXFLAGS_NO_UNUSED_FUNCTION "-Wno-unused-function" GMXC_CXXFLAGS)
        endif()
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
    endif()

    # Apple bastardized version of Clang
    if(${CMAKE_C_COMPILER_ID} MATCHES "AppleClang")
        if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 11.0)
            # Mac OS Catalina ships with a horribly broken compiler (version 11.0.0.11000033)
            # that checks stack alignment by default, but their own C library
            # does not align the stack properly. Embarrassing, Apple...
            GMX_TEST_CFLAG(CFLAG_NO_STACK_CHECK "-fno-stack-check" GMXC_CFLAGS)
        endif()
    endif()

    if(${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")
        if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 11.0)
            # Mac OS Catalina ships with a horribly broken compiler (version 11.0.0.11000033)
            # that checks stack alignment by default, but their own C library
            # does not align the stack properly. Embarrassing, Apple...
            GMX_TEST_CXXFLAG(CXXFLAG_NO_STACK_CHECK "-fno-stack-check" GMXC_CXXFLAGS)
        endif()
    endif()

    # Apple bastardized version of Clang
    if(${CMAKE_C_COMPILER_ID} MATCHES "AppleClang")
        if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 11.0)
            # Mac OS Catalina ships with a horribly broken compiler (version 11.0.0.11000033)
            # that checks stack alignment by default, but their own C library
            # does not align the stack properly. Embarrassing, Apple...
            GMX_TEST_CFLAG(CFLAG_NO_STACK_CHECK "-fno-stack-check" GMXC_CFLAGS)
        endif()
    endif()

    if(${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")
        if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 11.0)
            # Mac OS Catalina ships with a horribly broken compiler (version 11.0.0.11000033)
            # that checks stack alignment by default, but their own C library
            # does not align the stack properly. Embarrassing, Apple...
            GMX_TEST_CXXFLAG(CXXFLAG_NO_STACK_CHECK "-fno-stack-check" GMXC_CXXFLAGS)
        endif()
    endif()

    # Fujitsu compilers on PrimeHPC/Sparc64
    if(${CMAKE_C_COMPILER_ID} MATCHES Fujitsu OR
       (${CMAKE_C_COMPILER_ID} MATCHES unknown AND ${CMAKE_C_COMPILER} MATCHES ^fcc))
        GMX_TEST_CFLAG(CFLAG_GNUCOMPAT "-Xg;-w" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAG_OPT "-Kfast,reduction,swp,simd=2,uxsimd,fsimple;-x100" GMXC_CFLAGS)
    endif()

    if(${CMAKE_CXX_COMPILER_ID} MATCHES Fujitsu OR
       (${CMAKE_CXX_COMPILER_ID} MATCHES unknown AND ${CMAKE_CXX_COMPILER} MATCHES ^FCC))
        GMX_TEST_CXXFLAG(CXXFLAG_GNUCOMPAT "-Xg;-w" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAG_OPT "-Kfast,reduction,swp,simd=2,uxsimd,fsimple;-x100" GMXC_CXXFLAGS)
    endif()

endmacro()
