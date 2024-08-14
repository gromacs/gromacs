#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2009- The GROMACS Authors
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

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

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
    set(CFLAGS
            ${SIMD_C_FLAGS}
            ${EXTRA_C_FLAGS}
            ${GMXC_CFLAGS}
         PARENT_SCOPE)

    set(CXXFLAGS
            ${SIMD_CXX_FLAGS}
            ${EXTRA_CXX_FLAGS}
            ${GMXC_CXXFLAGS}
        PARENT_SCOPE)
endfunction()

# Implementation function to add compiler flags expected for all
# GROMACS build configurations, and those expected for the current
# CMake build type (e.g. Release) to TARGET. Other GROMACS CMake code
# is expected to use either gmx_target_compile_options(name_of_target)
# or gmx_device_target_compile_options(name_of_variable) because CUDA
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
    # TODO: Restrict the scope of MPI dependence.
    # Targets that actually need MPI headers and build tool flags should
    # manage their own `target_link_libraries` locally. Such a change is beyond
    # the scope of the bug fix for #4678.
    if (GMX_LIB_MPI AND TARGET ${TARGET})
        target_link_libraries(
            ${TARGET} PRIVATE
            $<$<LINK_LANGUAGE:CXX>:MPI::MPI_CXX>
            # We don't know whether we have sought the MPI::C component at all, or at least
            # by the time we process these lines.
            $<$<AND:$<TARGET_EXISTS:MPI::MPI_C>,$<LINK_LANGUAGE:C>>:MPI::MPI_C>
            $<$<LINK_LANGUAGE:CUDA>:MPI::MPI_CXX>
        )
    endif ()
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
# This function is also used for wrapping the options passed to the HIP
# compiler
function(gmx_device_target_compile_options VARIABLE_NAME)
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

function(gmx_target_interface_warning_suppression TARGET WARNING_FLAG VARNAME)
    check_cxx_compiler_flag(${WARNING_FLAG} ${VARNAME})
    if(${VARNAME})
        target_compile_options(${TARGET} INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${WARNING_FLAG}>)
        if (GMX_GPU_HIP)
            target_compile_options(${TARGET} INTERFACE $<$<COMPILE_LANGUAGE:HIP>:${WARNING_FLAG}>)
        endif()
    endif()
endfunction()

# This is the actual exported function to be called
macro (gmx_c_flags)

    include(CheckCCompilerFlag)
    include(CheckCXXCompilerFlag)

    # gcc
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        #flags are added in reverse order and -Wno* need to appear after -Wall
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall;-Wno-unused;-Wunused-value;-Wunused-parameter" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wextra;-Wno-sign-compare;-Wpointer-arith" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_UNDEF "-Wundef" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CFLAGS_RELEASE_ONLY)
            if(CYGWIN)
                GMX_TEST_CFLAG(CFLAGS_WARN_SUBSCRIPT "-Wno-char-subscripts" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_STRINGOP_TRUNCATION "-Werror=stringop-truncation" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN_NO_MISSING_FIELD_INITIALIZERS "-Wno-missing-field-initializers" GMXC_CFLAGS)
        # new in gcc 4.5
        GMX_TEST_CFLAG(CFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_COPT "-funroll-all-loops"
                       GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_NOINLINE "-fno-inline" GMXC_CFLAGS_DEBUG)
    endif()
    # g++
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            # Problematic with CUDA
            # GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EFFCXX "-Wnon-virtual-dtor" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra;-Wpointer-arith;-Wmissing-declarations" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_UNDEF "-Wundef" GMXC_CXXFLAGS)
            GMX_TEST_CFLAG(CXXFLAGS_WARN_REL "-Wno-array-bounds" GMXC_CXXFLAGS_RELEASE_ONLY)
            GMX_TEST_CXXFLAG(CXXFLAGS_STRINGOP_TRUNCATION "-Wstringop-truncation" GMXC_CXXFLAGS)
            if (CXXFLAGS_STRINGOP_TRUNCATION)
                set(CXXFLAGS_NO_STRINGOP_TRUNCATION "-Wno-stringop-truncation")
            endif()
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_MISSING_FIELD_INITIALIZERS "-Wno-missing-field-initializers" GMXC_CXXFLAGS)
        # new in gcc 4.5
        GMX_TEST_CXXFLAG(CXXFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_COPT "-funroll-all-loops"
                         GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_NOINLINE "-fno-inline" GMXC_CXXFLAGS_DEBUG)
    endif()
    # NVHPC Compiler supported from version 24.1.
    if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
        GMX_TEST_CXXFLAG(CXXFLAGS_NVCXX_FAST "-fast" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_NVCXX_NOFLUSHZ "-Mnoflushz" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_NVCXX_NODAZ "-Mnodaz" GMXC_CXXFLAGS)
        GMX_TEST_CFLAG(CFLAGS_NVC_FAST "-fast" GMXC_CFLAGS)
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
        GMX_TEST_CXXFLAG(CXXFLAGS_LANG "/permissive-" GMXC_CXXFLAGS)
        # Set source and execution character sets to UTF-8
        GMX_TEST_CXXFLAG(CXXFLAGS_LANG "/utf-8" GMXC_CXXFLAGS)
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang" OR CMAKE_C_COMPILER_ID MATCHES "IntelLLVM")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        if (GMX_COMPILER_WARNINGS)
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall;-Wno-unused;-Wunused-value;-Wunused-parameter" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wpointer-arith" GMXC_CFLAGS_EXTRA)
        endif()
        if (CMAKE_BUILD_TYPE MATCHES "Debug")
            GMX_TEST_CFLAG(CFLAGS_NO_DEBUG_DISABLES_OPTIMIZATION "-Wno-debug-disables-optimization" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN_NO_MISSING_FIELD_INITIALIZERS "-Wno-missing-field-initializers" GMXC_CFLAGS)
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        if (GMX_COMPILER_WARNINGS)
            # If used, -Wall should precede other options that silence warnings it enables
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.0") #LLVM BUG #21629
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_BRACES "-Wno-missing-braces" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra;-Wpointer-arith;-Wmissing-prototypes" GMXC_CXXFLAGS)
            if(NOT CMAKE_CXX_FLAGS MATCHES "-Wno-deprecated")
                GMX_TEST_CXXFLAG(CXXFLAGS_DEPRECATED "-Wdeprecated" GMXC_CXXFLAGS)
            endif()

            if (APPLE)
                # macOS Ventura deprecated `sprintf` in favor of `snprintf`.
                # This workaround suppresses the deprecation warnings.
                GMX_TEST_CXXFLAG(CXXFLAGS_NO_DEPRECATED_DECLARATIONS "-Wno-deprecated-declarations" GMXC_CXXFLAGS)
                # This warning is only useful for cross-compiling
                GMX_TEST_CXXFLAG(CXXFLAGS_NO_POISON_SYSTEM_DIRECTORIES "-Wno-poison-system-directories" GMXC_CXXFLAGS)
            endif()

            # Functions placed in headers for inlining are not always
            # used in every translation unit that includes the files,
            # so we must disable the warning that there are such
            # functions that are unused.
            GMX_TEST_CXXFLAG(CXXFLAGS_NO_UNUSED_FUNCTION "-Wno-unused-function" GMXC_CXXFLAGS)
        endif()
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_RESERVED_IDENTIFIER "-Wno-reserved-identifier" GMXC_CXXFLAGS) # LLVM BUG #50644
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_MISSING_FIELD_INITIALIZERS "-Wno-missing-field-initializers" GMXC_CXXFLAGS)
        if(GMX_GPU)
            string(TOUPPER "${GMX_GPU}" _gmx_gpu_uppercase)
            if(${_gmx_gpu_uppercase} STREQUAL "CUDA")
                # disable the missing prototypes warning as there are some tests which builds with nvcc and clang as host compiler
                # still have this issue
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_MISSING_PROTOTYPES "-Wno-missing-prototypes" GMXC_CXXFLAGS)
                GMX_TEST_CXXFLAG(CXXFLAGS_WARN_NO_CONSTANT_LOGICAL_OPERAND "-Wno-constant-logical-operand" GMXC_CXXFLAGS)
            endif()
        endif()

        if (CMAKE_BUILD_TYPE MATCHES "Debug")
            GMX_TEST_CXXFLAG(CXXFLAGS_NO_DEBUG_DISABLES_OPTIMIZATION "-Wno-debug-disables-optimization" GMXC_CXXFLAGS)
        endif()
    endif()

    # Apple bastardized version of Clang
    if(${CMAKE_C_COMPILER_ID} MATCHES "AppleClang")
        if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 11.0)
            # macOS Catalina ships with a horribly broken compiler (version 11.0.0.11000033)
            # that checks stack alignment by default, but their own C library
            # does not align the stack properly. Embarrassing, Apple...
            GMX_TEST_CFLAG(CFLAG_NO_STACK_CHECK "-fno-stack-check" GMXC_CFLAGS)
        endif()
    endif()

    if(${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")
        if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 11.0)
            # macOS Catalina ships with a horribly broken compiler (version 11.0.0.11000033)
            # that checks stack alignment by default, but their own C library
            # does not align the stack properly. Embarrassing, Apple...
            GMX_TEST_CXXFLAG(CXXFLAG_NO_STACK_CHECK "-fno-stack-check" GMXC_CXXFLAGS)
        endif()
    endif()

endmacro()

# Make sure we generate warnings (and hopefully fix) "everything"
# reported by compilers that support that (ie recent clang and its
# derivatives).
function(gmx_warn_on_everything target)

    # If the compiler suports warning on "everything" then we'll turn
    # it on. Note that all warnings become errors for developer
    # builds, but not for user builds.
    gmx_target_warning_suppression(${target} "-Weverything" HAS_WARNING_EVERYTHING)

    if (NOT HAS_WARNING_EVERYTHING)
        # There's no need to suppress aspects of "-Weverything" if
        # that warning is not supported.
        return()
    endif()

    # We don't actually fix everything, so list the exceptions that we
    # choose to make. We may be able to eliminate some of these over
    # time.
    #
    # We check whether the flag is accepted first, so that we suppress
    # such warnings also with compilers that don't directly identify
    # as e.g. clang despite being based on it (e.g. most vendor
    # compilers), and also don't fail to compile GROMACS when future
    # versions of any such compiler changes how the warnings
    # look/work.

    # We have no intention of C++98 compatibility
    gmx_target_warning_suppression(${target} "-Wno-c++98-compat" HAS_WARNING_NO_CPLUSPLUS98_COMPAT)
    gmx_target_warning_suppression(${target} "-Wno-c++98-compat-pedantic" HAS_WARNING_NO_CPLUSPLUS98_COMPAT_PEDANTIC)

    # We require newer C++ than version 11, so ignore c++11 specific warnings
    gmx_target_warning_suppression(${target} "-Wno-return-std-move-in-c++11" HAS_WARNING_NO_RETURN_STD_MOVE_IN_CPLUSPLUS11)

    # Don't warn for use of OpenMP pragmas in no-omp build
    gmx_target_warning_suppression(${target} "-Wno-source-uses-openmp" HAS_WARNING_NO_SOURCE_USED_OPENMP)

    # Allowed in attributes (compilers are required to ignore unknown attributes)
    gmx_target_warning_suppression(${target} "-Wno-c++17-extensions" HAS_WARNING_NO_CPLUSPLUS17_EXTENSIONS)

    # Custom Doxygen commands are used
    gmx_target_warning_suppression(${target} "-Wno-documentation-unknown-command" HAS_WARNING_NO_DOCUMENTATION_UNKNOWN_COMMAND)

    # We need to use default labels in switch statements, because GCC gives
    # maybe-uninitialized without default label and checks for illegal enum values.
    gmx_target_warning_suppression(${target} "-Wno-covered-switch-default" HAS_WARNING_NO_COVERED_SWITCH_DEFAULT)

    # Default statement for enum is OK.
    # It's OK to not have branches for Count members of enum classes
    gmx_target_warning_suppression(${target} "-Wno-switch-enum" HAS_WARNING_NO_SWITCH_ENUM)

    # Sometimes it's ok to not have a default: switch statement
    gmx_target_warning_suppression(${target} "-Wno-switch-default" HAS_WARNING_NO_SWITCH_DEFAULT)

    # We need to use macros like
    # GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR and
    # CLANG_DIAGNOSTIC_IGNORE. Those will look strange if they don't
    # have a semicolon after them, and might confuse tools like IDEs
    # also.
    gmx_target_warning_suppression(${target} "-Wno-extra-semi-stmt" HAS_WARNING_NO_EXTRA_SEMI_STMT)

    # We intend to use fully inline classes with virtual methods
    gmx_target_warning_suppression(${target} "-Wno-weak-vtables" HAS_WARNING_NO_WEAK_VTABLES)

    # We intend to use constructor arguments that shadow member variables
    gmx_target_warning_suppression(${target} "-Wno-shadow" HAS_WARNING_NO_SHADOW)

    # Padding of structs is routine, we don't need to hear about it
    gmx_target_warning_suppression(${target} "-Wno-padded" HAS_WARNING_NO_PADDED)

    # Our uses of double underscores in macro names are OK
    gmx_target_warning_suppression(${target} "-Wno-reserved-id-macro" HAS_WARNING_NO_RESERVED_ID_MACRO)

    # Implicit conversion of float to double is fine
    gmx_target_warning_suppression(${target} "-Wno-double-promotion" HAS_WARNING_NO_DOUBLE_PROMOTION)

    # No resources in static variables need exit-time destructors
    gmx_target_warning_suppression(${target} "-Wno-exit-time-destructors" HAS_WARNING_NO_EXIT_TIME_DESTRUCTORS)

    # Global constructors are not needed
    gmx_target_warning_suppression(${target} "-Wno-global-constructors" HAS_WARNING_NO_GLOBAL_CONSTRUCTORS)

    # False positives are emitted
    gmx_target_warning_suppression(${target} "-Wno-documentation" HAS_WARNING_NO_DOCUMENTATION)

    # We intend to use format strings that we construct, even though that is a security risk
    gmx_target_warning_suppression(${target} "-Wno-format-nonliteral" HAS_WARNING_NO_FORMAT_NONLITERAL)

    # We do a lot of conditional compilation that sometimes uses a symbol and sometimes does not
    gmx_target_warning_suppression(${target} "-Wno-used-but-marked-unused" HAS_WARNING_NO_USED_BUT_MARKED_UNUSED)

    # It's only risky to compare floats for equality when they are the
    # result of computation.  Unfortunately it's hard to tell the
    # difference and there's no good way to suppress this on a
    # case-by-base basis.
    gmx_target_warning_suppression(${target} "-Wno-float-equal" HAS_WARNING_NO_FLOAT_EQUAL)

    if(GMX_GPU_SYCL)
        # Intel oneAPI compiler warns about compiling kernels with non-32 subgroups for CUDA devices.
        # But we can not avoid compiling kernels for other sub-group sizes unless we have
        # sycl_ext_oneapi_kernel_properties and/or sycl::any_device_has (https://github.com/intel/llvm/issues/5562).
        # Even then, we might want to build for multiple devices. So, we silence the warning.
        gmx_target_warning_suppression(${target} "-Wno-cuda-compat" HAS_WARNING_NO_CUDA_COMPAT)
    endif()

    #
    # Exceptions we should consider fixing
    #

    # Much code in gmxana uses complex logic that may or may not be valid
    gmx_target_warning_suppression(${target} "-Wno-conditional-uninitialized" HAS_WARNING_CONDITIONAL_UNINITIALIZED)

    # We have many places implicit conversions still occur, most of which need fixing
    gmx_target_warning_suppression(${target} "-Wno-conversion" HAS_WARNING_NO_CONVERSION)

    # We use the Linux signal handlers in the intended way, but it triggers this warning.
    # It would be better to localize this exception.
    gmx_target_warning_suppression(${target} "-Wno-disabled-macro-expansion" HAS_WARNING_NO_DISABLED_MACRO_EXPANSION)

    # The NBNXM simd kernels define lots of macros that are not used
    # It would be better to localize this exception.
    gmx_target_warning_suppression(${target} "-Wno-unused-macros" HAS_WARNING_NO_UNUSED_MACROS)

    # The C++ Buffer hardening warning in Clang is still experimental as of January 2023.
    # After it is more mature it would better to localize this exception.
    gmx_target_warning_suppression(${target} "-Wno-unsafe-buffer-usage" HAS_WARNING_NO_UNSAFE_BUFFER_USAGE)

endfunction()
