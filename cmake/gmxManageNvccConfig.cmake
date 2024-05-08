#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
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
# away from users installing GROMACS. See Issue #1982.
function(work_around_glibc_2_23)
    try_compile(IS_GLIBC_2_23_OR_HIGHER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/TestGlibcVersion.cpp)
    if(IS_GLIBC_2_23_OR_HIGHER)
        message(STATUS "Adding work-around for issue compiling CUDA code with glibc 2.23 string.h")
        list(APPEND CUDA_HOST_COMPILER_OPTIONS "-D_FORCE_INLINES")
        set(CUDA_HOST_COMPILER_OPTIONS ${CUDA_HOST_COMPILER_OPTIONS} PARENT_SCOPE)
    endif()
endfunction()

gmx_check_if_changed(CUDA_HOST_COMPILER_CHANGED CMAKE_CUDA_HOST_COMPILER)

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

    mark_as_advanced(CMAKE_CUDA_HOST_COMPILER  CUDA_HOST_COMPILER_OPTIONS)
endif()

# We would like to be helpful and reject the host compiler with a
# clear error message at configure time, rather than let nvcc
# later reject the host compiler as not supported when the first
# CUDA source file is built. We've implemented that for current
# nvcc running on Unix-like systems, but e.g. changes to nvcc
# will further affect the limited portability of this checking
# code. Set the CMake variable GMX_NVCC_WORKS on if you want to
# bypass this check.
if((_cuda_nvcc_executable_or_flags_changed OR CUDA_HOST_COMPILER_CHANGED OR NOT GMX_NVCC_WORKS) AND NOT WIN32)
    message(STATUS "Check for working NVCC/C++ compiler combination with nvcc '${CUDAToolkit_NVCC_EXECUTABLE}'")

    execute_process(COMMAND ${CUDAToolkit_NVCC_EXECUTABLE} --compiler-bindir=${CMAKE_CUDA_HOST_COMPILER} -c ${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS_${_build_type}} ${CMAKE_SOURCE_DIR}/cmake/TestCUDA.cu
        RESULT_VARIABLE _cuda_test_res
        OUTPUT_VARIABLE _cuda_test_out
        ERROR_VARIABLE  _cuda_test_err
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(${_cuda_test_res})
        message(STATUS "Check for working NVCC/C compiler combination - broken")
        message(STATUS "${CUDAToolkit_NVCC_EXECUTABLE} standard output: '${_cuda_test_out}'")
        message(STATUS "${CUDAToolkit_NVCC_EXECUTABLE} standard error:  '${_cuda_test_err}'")
        if(${_cuda_test_err} MATCHES "nsupported")
            message(FATAL_ERROR "NVCC/C++ compiler combination does not seem to be supported. CUDA frequently does not support the latest versions of the host compiler, so you might want to try an earlier C++ compiler version and make sure your CUDA compiler and driver are as recent as possible. Set the GMX_NVCC_WORKS CMake cache variable to bypass this check if you know what you are doing.")
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 11.2 AND CUDAToolkit_VERSION VERSION_GREATER 11.4 AND CUDAToolkit_VERSION VERSION_LESS 11.7) # Issue #4574, #4641
            # Above, we should be checking for VERSION_LESS 11.6.2, but CUDA_VERSION is only "major.minor"
            message(FATAL_ERROR "CUDA versions 11.5-11.6.1 are known to be incompatible with some GCC 11.x. Use a different GCC or update your CUDA installation to at least CUDA 11.6.2")
        else()
            message(FATAL_ERROR "CUDA compiler does not seem to be functional or is not compatible with the host compiler. Set the GMX_NVCC_WORKS CMake cache variable to bypass this check if you know what you are doing.")
        endif()
    elseif(NOT GMX_CUDA_TEST_COMPILER_QUIETLY)
        message(STATUS "Check for working NVCC/C++ compiler combination - works")
        set(GMX_NVCC_WORKS TRUE CACHE INTERNAL "Nvcc can compile a trivial test program")
    endif()
endif() # GMX_CHECK_NVCC

# Tests a single flag to use with nvcc.
#
# If the flags are accepted, they are appended to the variable named
# in the first argument. The cache variable named in the second
# argument is used to avoid rerunning the check in future invocations
# of cmake. The list of flags to check follows these two required
# arguments.
#
# Note that a space-separated string of flags, or a flag-value pair
# separated by spaces will not work. Use the single-argument forms
# accepted by nvcc, like "--arg=value".
#
# As this code is not yet tested on Windows, it always accepts the
# flags in that case.
function(gmx_add_nvcc_flag_if_supported _output_variable_name_to_append_to _flags_cache_variable_name)
    # If the check has already been run, do not re-run it
    if (NOT ${_flags_cache_variable_name} AND NOT WIN32)
        message(STATUS "Checking if nvcc accepts flags ${ARGN}")
        execute_process(
            COMMAND ${CUDAToolkit_NVCC_EXECUTABLE} ${ARGN} -ccbin ${CMAKE_CUDA_HOST_COMPILER} "${CMAKE_SOURCE_DIR}/cmake/TestCUDA.cu"
            RESULT_VARIABLE _cuda_success
            OUTPUT_QUIET
            ERROR_QUIET
            )
        # Convert the success value to a boolean and report status
        if (_cuda_success EQUAL 0)
            set(_cache_variable_value TRUE)
            message(STATUS "Checking if nvcc accepts flags ${ARGN} - Success")
        else()
            if(NOT(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11))
              set(CCBIN "--compiler-bindir=${CMAKE_CUDA_HOST_COMPILER}")
            endif()
            execute_process(
                COMMAND ${CUDAToolkit_NVCC_EXECUTABLE} ${ARGN} ${CCBIN} "${CMAKE_SOURCE_DIR}/cmake/TestCUDA.cu"
                RESULT_VARIABLE _cuda_success
                OUTPUT_QUIET
                ERROR_QUIET
                )
            # Convert the success value to a boolean and report status
            if (_cuda_success EQUAL 0)
                set(_cache_variable_value TRUE)
                message(STATUS "Checking if nvcc accepts flags ${ARGN} - Yes")
            else()
                set(_cache_variable_value FALSE)
                message(STATUS "Checking if nvcc accepts flags ${ARGN} - No")
            endif()
        endif()
        set(${_flags_cache_variable_name} ${_cache_variable_value} CACHE INTERNAL "Whether NVCC supports flag(s) ${ARGN}")
    endif()
    # Append the flags to the output variable if they have been tested to work
    if (${_flags_cache_variable_name} OR WIN32)
        list(APPEND ${_output_variable_name_to_append_to} ${ARGN})
        set(${_output_variable_name_to_append_to} ${${_output_variable_name_to_append_to}} PARENT_SCOPE)
    endif()
endfunction()

# If any of these manual override variables for target CUDA GPU architectures
# or virtual architecture is set, parse the values and assemble the nvcc
# command line for these. Otherwise use our defaults.
# Note that the manual override variables require a semicolon separating
# architecture codes.
set(GMX_CUDA_NVCC_GENCODE_FLAGS)
if (GMX_CUDA_TARGET_SM OR GMX_CUDA_TARGET_COMPUTE)
    set(_target_sm_list ${GMX_CUDA_TARGET_SM})
    foreach(_target ${_target_sm_list})
        gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_${_target} "--generate-code=arch=compute_${_target},code=sm_${_target}")
        if (NOT NVCC_HAS_GENCODE_COMPUTE_AND_SM_${_target} AND NOT WIN32)
            message(FATAL_ERROR "Your choice of ${_target} in GMX_CUDA_TARGET_SM was not accepted by nvcc, please choose a target that it accepts")
        endif()
    endforeach()
    set(_target_compute_list ${GMX_CUDA_TARGET_COMPUTE})
    foreach(_target ${_target_compute_list})
        gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_${_target} --generate-code=arch=compute_${_target},code=compute_${_target})
        if (NOT NVCC_HAS_GENCODE_COMPUTE_${_target} AND NOT WIN32)
            message(FATAL_ERROR "Your choice of ${_target} in GMX_CUDA_TARGET_COMPUTE was not accepted by nvcc, please choose a target that it accepts")
        endif()
    endforeach()
else()
    # TODO: Note that this SM arch populating code can be replaced with `CUDA_ARCHITECTURES=all` once we upgrade cmake to >= 3.23

    # Set the CUDA GPU architectures to compile for:
    # - with CUDA >=11.0        CC 8.0 is supported
    #     => compile sm_35, sm_37, sm_50, sm_52, sm_60, sm_61, sm_70, sm_75, sm_80 SASS, and compute_35, compute_80 PTX

    # First add flags that trigger SASS (binary) code generation for physical arch
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_35 --generate-code=arch=compute_35,code=sm_35)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_37 --generate-code=arch=compute_37,code=sm_37)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_50 --generate-code=arch=compute_50,code=sm_50)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_52 --generate-code=arch=compute_52,code=sm_52)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_60 --generate-code=arch=compute_60,code=sm_60)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_61 --generate-code=arch=compute_61,code=sm_61)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_70 --generate-code=arch=compute_70,code=sm_70)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_75 --generate-code=arch=compute_75,code=sm_75)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_80 --generate-code=arch=compute_80,code=sm_80)
    # Don't attempt to add newest architectures with old GNU compiler, to avoid issues in CI
    # related to being unable to test which flags are supported
    if (NOT (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8))
        gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_86 --generate-code=arch=compute_86,code=sm_86)
        gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_89 --generate-code=arch=compute_89,code=sm_89)
        gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_AND_SM_90 --generate-code=arch=compute_90,code=sm_90)
    endif()
    # Requesting sm or compute 35, 37, or 50 triggers deprecation messages with
    # nvcc 11.0, which we need to suppress for use in CI
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_FLAGS NVCC_HAS_WARNING_NO_DEPRECATED_GPU_TARGETS -Wno-deprecated-gpu-targets)

    # Next add flags that trigger PTX code generation for the
    # newest supported virtual arch that's useful to JIT to future architectures
    # as well as an older one suitable for JIT-ing to any rare intermediate arch
    # (like that of Jetson / Drive PX devices)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_53 --generate-code=arch=compute_53,code=sm_53)
    gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_GENCODE_FLAGS NVCC_HAS_GENCODE_COMPUTE_80 --generate-code=arch=compute_80,code=sm_80)
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
gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_FLAGS NVCC_HAS_USE_FAST_MATH -use_fast_math)
# Add warnings
gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_FLAGS NVCC_HAS_PTXAS_WARN_DOUBLE_USAGE -Xptxas=-warn-double-usage)
gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_FLAGS NVCC_HAS_PTXAS_WERROR -Xptxas=-Werror)

# strip gencode/arch from GMX_CUDA_NVCC_GENCODE_FLAGS only keep the arch numbers
string(REGEX REPLACE "([A-Za-z=_-])" ""  GMX_CUDA_NVCC_GENCODE_FLAGS "${GMX_CUDA_NVCC_GENCODE_FLAGS}")
string(REGEX REPLACE "([0-9]+,)" ""  GMX_CUDA_NVCC_GENCODE_FLAGS "${GMX_CUDA_NVCC_GENCODE_FLAGS}")

# Set the openmp host compilation flag if it has not been set automatically.
# Some other compilation flags, mostly warnings are not automatically
# propagated, but we do not set them all manually. See issue #4757 and the
# MR discussion in !3780 for more information.
if(GMX_OPENMP AND NOT "${OpenMP_CXX_FLAGS}" STREQUAL "")
    if(NOT ${OpenMP_CXX_FLAGS} IN_LIST GMX_CUDA_NVCC_FLAGS)
        list(APPEND CUDA_HOST_COMPILER_OPTIONS -Xcompiler=${OpenMP_CXX_FLAGS})
    endif()
endif()

# pass the CMAKE_CXX_FLAGS to nvcc as host compiler options
if(NOT "${CMAKE_CXX_FLAGS}" STREQUAL "")
    if(NOT ${CMAKE_CXX_FLAGS} IN_LIST GMX_CUDA_NVCC_FLAGS)
        # convert the string to list
        string(REPLACE " " ";" CMAKE_CXX_FLAGS_LIST_ "${CMAKE_CXX_FLAGS}")
        # Loop over each item in the list
        foreach(item IN LISTS CMAKE_CXX_FLAGS_LIST_)
            list(APPEND CUDA_HOST_COMPILER_OPTIONS -Xcompiler="${item}")
        endforeach()
    endif()
endif()

# pass the GMXC_CXXFLAGS to nvcc as host compiler options.
# This is used to suppress warnings for Clang for now.
# Ideally should be enabled for all host compilers once
# all the warnings are resolved in the test cases which are
# compiled with nvcc.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(NOT "${GMXC_CXXFLAGS}" STREQUAL "")
        if(NOT GMXC_CXXFLAGS IN_LIST GMX_CUDA_NVCC_FLAGS)
            # convert the string to list
            string(REPLACE " " ";" GMXC_CXXFLAGS_LIST_ "${GMXC_CXXFLAGS}")
            # Loop over each item in the list
            foreach(item IN LISTS GMXC_CXXFLAGS_LIST_)
                list(APPEND CUDA_HOST_COMPILER_OPTIONS -Xcompiler="${item}")
            endforeach()
        endif()
    endif()
endif()

# assemble the CUDA host compiler flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${CUDA_HOST_COMPILER_OPTIONS}")


# Disable cudafe warnings with nvc++ as a host compiler - warning #177-D
gmx_add_nvcc_flag_if_supported(GMX_CUDA_NVCC_FLAGS NVCC_HAS_DIAG_SUPPRESS_177 -diag-suppress=177)

string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type)
gmx_check_if_changed(_cuda_nvcc_executable_or_flags_changed CUDAToolkit_NVCC_EXECUTABLE CUDA_NVCC_FLAGS CUDA_NVCC_FLAGS_${_build_type})
