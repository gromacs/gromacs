#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
# Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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

# Manage CUDA nvcc compilation configuration, try to be smart to ease the users'
# pain as much as possible:
# - use the CUDA_HOST_COMPILER if defined by the user, otherwise
# - check if nvcc works with CUDA_HOST_COMPILER and the generated nvcc and C++ flags
#
# - (advanced) variables set:
#   * CUDA_HOST_COMPILER_OPTIONS    - the full host-compiler related option list passed to nvcc
#
# Note that from CMake 2.8.10 FindCUDA defines CUDA_HOST_COMPILER internally,
# so we won't set it ourselves, but hope that the module does a good job.

# glibc 2.23 changed string.h in a way that breaks CUDA compilation in
# many projects, but which has a trivial workaround. It would be nicer
# to compile with nvcc and see that the workaround is necessary and
# effective, but it is unclear how to do that. Also, grepping in the
# glibc source shows that _FORCE_INLINES is only used in this string.h
# feature and performance of memcpy variants is unimportant for CUDA
# code in GROMACS. So this workaround is good enough to keep problems
# away from users installing GROMACS. See Issue #1942.
function(work_around_glibc_2_23)
    try_compile(IS_GLIBC_2_23_OR_HIGHER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/TestGlibcVersion.cpp)
    if(IS_GLIBC_2_23_OR_HIGHER)
        message(STATUS "Adding work-around for issue compiling CUDA code with glibc 2.23 string.h")
        list(APPEND CUDA_HOST_COMPILER_OPTIONS "-D_FORCE_INLINES")
        set(CUDA_HOST_COMPILER_OPTIONS ${CUDA_HOST_COMPILER_OPTIONS} PARENT_SCOPE)
    endif()
endfunction()

gmx_check_if_changed(CUDA_HOST_COMPILER_CHANGED CUDA_HOST_COMPILER)

# set up host compiler and its options
if(CUDA_HOST_COMPILER_CHANGED)
    set(CUDA_HOST_COMPILER_OPTIONS "")

    if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "GNU")
        # Some versions of gcc-4.8 and gcc-4.9 have produced errors
        # (in particular on OS X) if we do not use
        # -D__STRICT_ANSI__. It is harmless, so we might as well add
        # it for all versions.
        list(APPEND CUDA_HOST_COMPILER_OPTIONS "-D__STRICT_ANSI__")
    endif()

    work_around_glibc_2_23()

    set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS}"
        CACHE STRING "Options for nvcc host compiler (do not edit!).")

    mark_as_advanced(CUDA_HOST_COMPILER CUDA_HOST_COMPILER_OPTIONS)
endif()

# If any of these manual override variables for target CUDA GPU architectures
# or virtual architecture is set, parse the values and assemble the nvcc
# command line for these. Otherwise use our defaults.
# Note that the manual override variables require a semicolon separated
# architectures codes.
if (GMX_CUDA_TARGET_SM OR GMX_CUDA_TARGET_COMPUTE)
    set(GMX_CUDA_NVCC_GENCODE_FLAGS)
    set(_target_sm_list ${GMX_CUDA_TARGET_SM})
    foreach(_target ${_target_sm_list})
        list(APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_${_target},code=sm_${_target}")
    endforeach()
    set(_target_compute_list ${GMX_CUDA_TARGET_COMPUTE})
    foreach(_target ${_target_compute_list})
        list(APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_${_target},code=compute_${_target}")
    endforeach()
else()
    # Set the CUDA GPU architectures to compile for:
    # - with CUDA >=9.0         CC 7.0 is supported and CC 2.0 is no longer supported
    #     => compile sm_30, sm_35, sm_37, sm_50, sm_52, sm_60, sm_61, sm_70 SASS, and compute_70 PTX
    # - with CUDA >=10.0        CC 7.5 is supported
    #     => compile sm_30, sm_35, sm_37, sm_50, sm_52, sm_60, sm_61, sm_70, sm_75 SASS, and compute_75 PTX

    # First add flags that trigger SASS (binary) code generation for physical arch
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_30,code=sm_30")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=sm_35")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_37,code=sm_37")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=sm_50")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_52,code=sm_52")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_60,code=sm_60")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_61,code=sm_61")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_70,code=sm_70")

    # Next add flags that trigger PTX code generation for the newest supported virtual arch
    # that's useful to JIT to future architectures
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=compute_35")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=compute_50")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_52,code=compute_52")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_60,code=compute_60")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_61,code=compute_61")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_70,code=compute_70")
    if(NOT CUDA_VERSION VERSION_LESS "10.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_75,code=compute_75")
    endif()
endif()

if (GMX_CUDA_TARGET_SM)
    set_property(CACHE GMX_CUDA_TARGET_SM PROPERTY HELPSTRING "List of CUDA GPU architecture codes to compile for (without the sm_ prefix)")
    set_property(CACHE GMX_CUDA_TARGET_SM PROPERTY TYPE STRING)
endif()
if (GMX_CUDA_TARGET_COMPUTE)
    set_property(CACHE GMX_CUDA_TARGET_COMPUTE PROPERTY HELPSTRING "List of CUDA virtual architecture codes to compile for (without the compute_ prefix)")
    set_property(CACHE GMX_CUDA_TARGET_COMPUTE PROPERTY TYPE STRING)
endif()

# assemble the CUDA flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_GENCODE_FLAGS}")
list(APPEND GMX_CUDA_NVCC_FLAGS "-use_fast_math")

# assemble the CUDA host compiler flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${CUDA_HOST_COMPILER_OPTIONS}")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type)
gmx_check_if_changed(_cuda_nvcc_executable_or_flags_changed CUDA_NVCC_EXECUTABLE CUDA_NVCC_FLAGS CUDA_NVCC_FLAGS_${_build_type})

# We would like to be helpful and reject the host compiler with a
# clear error message at configure time, rather than let nvcc
# later reject the host compiler as not supported when the first
# CUDA source file is built. We've implemented that for current
# nvcc running on Unix-like systems, but e.g. changes to nvcc
# will further affect the limited portability of this checking
# code. Set the CMake variable GMX_NVCC_WORKS on if you want to
# bypass this check.
if((_cuda_nvcc_executable_or_flags_changed OR CUDA_HOST_COMPILER_CHANGED OR NOT GMX_NVCC_WORKS) AND NOT WIN32)
    message(STATUS "Check for working NVCC/C++ compiler combination with nvcc '${CUDA_NVCC_EXECUTABLE}'")
    execute_process(COMMAND ${CUDA_NVCC_EXECUTABLE} -ccbin ${CUDA_HOST_COMPILER} -c ${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS_${_build_type}} ${CMAKE_SOURCE_DIR}/cmake/TestCUDA.cu
        RESULT_VARIABLE _cuda_test_res
        OUTPUT_VARIABLE _cuda_test_out
        ERROR_VARIABLE  _cuda_test_err
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(${_cuda_test_res})
        message(STATUS "Check for working NVCC/C compiler combination - broken")
        message(STATUS "${CUDA_NVCC_EXECUTABLE} standard output: '${_cuda_test_out}'")
        message(STATUS "${CUDA_NVCC_EXECUTABLE} standard error:  '${_cuda_test_err}'")
        if(${_cuda_test_err} MATCHES "nsupported")
            message(FATAL_ERROR "NVCC/C++ compiler combination does not seem to be supported. CUDA frequently does not support the latest versions of the host compiler, so you might want to try an earlier C++ compiler version and make sure your CUDA compiler and driver are as recent as possible.")
        else()
            message(FATAL_ERROR "CUDA compiler does not seem to be functional.")
        endif()
    elseif(NOT GMX_CUDA_TEST_COMPILER_QUIETLY)
        message(STATUS "Check for working NVCC/C++ compiler combination - works")
        set(GMX_NVCC_WORKS TRUE CACHE INTERNAL "Nvcc can compile a trivial test program")
    endif()
endif() # GMX_CHECK_NVCC


# The flags are set as local variables which shadow the cache variables. The cache variables
# (can be set by the user) are appended. This is done in a macro to set the flags when all
# host compiler flags are already set.
macro(GMX_SET_CUDA_NVCC_FLAGS)
    set(CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_FLAGS};${CUDA_NVCC_FLAGS}")
endmacro()

# This helper function creates a temporary scope in which we can set
# the definitions, include directories and magic host-compiler-flag
# variables that have to be set in advance of calling
# cuda_add_library(). This is the only way cuda_add_library() can
# modify the command line used for host compilation. It is not
# possible to use the standard CMake mechanisms like
# target_compile_options() to add such things to targets after they
# are created.
function(gmx_cuda_add_library TARGET)
    add_definitions(-DHAVE_CONFIG_H)
    # Source files generated by NVCC can include gmxmpi.h, and so
    # need access to thread-MPI.
    include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
    # Source files can also contain topology related files and need access to
    # the remaining external headers
    include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/src/external)

    # Now add all the compilation options
    gmx_cuda_target_compile_options(CUDA_${TARGET}_CXXFLAGS)
    list(APPEND CMAKE_CXX_FLAGS ${CUDA_${TARGET}_CXXFLAGS})
    foreach(build_type ${build_types_with_explicit_flags})
        list(APPEND CMAKE_CXX_FLAGS_${build_type} ${CUDA_${TARGET}_CXXFLAGS_${build_type}})
    endforeach()

    cuda_add_library(${TARGET} ${ARGN})
endfunction()
