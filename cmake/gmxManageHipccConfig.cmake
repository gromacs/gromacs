#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2024- The GROMACS Authors
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

# Manage HIP clang compilation configuration, try to be smart to ease the users'
# pain as much as possible:

# Tests a single flag to use with hipcc.
#
# If the flags are accepted, they are appended to the variable named
# in the first argument. The cache variable named in the second
# argument is used to avoid rerunning the check in future invocations
# of cmake. The list of flags to check follows these two required
# arguments.
#
# Note that a space-separated string of flags, or a flag-value pair
# separated by spaces will not work. Use the single-argument forms
# accepted by hipcc, like "--arg=value".
#
# As this code is not yet tested on Windows, it always accepts the
# flags in that case.
#
# To avoid triggering device autodetection on systems that have no
# device present, we always add one architecture to the compilation
# arguments.
function(gmx_add_hipcc_flag_if_supported _output_variable_name_to_append_to _flags_cache_variable_name)
    # If the check has already been run, do not re-run it
    if (NOT DEFINED ${_flags_cache_variable_name} AND NOT WIN32)
        message(STATUS "Checking if hipcc accepts flags ${ARGN}")
        execute_process(
                COMMAND ${HIP_HIPCC_EXECUTABLE} ${ARGN} --offload-arch=gfx90a -Werror "${CMAKE_SOURCE_DIR}/cmake/TestHIP.cpp"
            RESULT_VARIABLE _hip_success
            OUTPUT_QUIET
            ERROR_QUIET
            )
        # Convert the success value to a boolean and report status
        if (_hip_success EQUAL 0)
            set(_cache_variable_value TRUE)
            message(STATUS "Checking if hipcc accepts flags ${ARGN} - Success")
        else()
            set(_cache_variable_value FALSE)
            message(STATUS "Checking if hipcc accepts flags ${ARGN} - No")
        endif()
        set(${_flags_cache_variable_name} ${_cache_variable_value} CACHE INTERNAL "Whether HIPCC supports flag(s) ${ARGN}")
    endif()
    # Append the flags to the output variable if they have been tested to work
    if (${_flags_cache_variable_name} OR WIN32)
        list(APPEND ${_output_variable_name_to_append_to} "${ARGN}")
        set(${_output_variable_name_to_append_to} ${${_output_variable_name_to_append_to}} PARENT_SCOPE)
    endif()
endfunction()

function(gmx_hip_check_single_flag _single_flag)
    string(REGEX REPLACE "=" "_" _flag_name_sanitized HIPCC_SUPPORTS_FLAG_${_single_flag})
    gmx_add_hipcc_flag_if_supported(GMX_HIP_HIPCC_FLAGS ${_flag_name_sanitized} ${_single_flag})
    if (NOT ${_flag_name_sanitized})
        message(STATUS "The version of the HIPCC compiler does not support ${_single_flag}")
    endif()
    set(GMX_HIP_HIPCC_FLAGS ${GMX_HIP_HIPCC_FLAGS} PARENT_SCOPE)
endfunction()

# iterate over a list of proposed target architectures and check that all of them are accepted by the compiler
# takes a list of target architectures
function(gmx_check_hip_architectures _target_architectures)
    string(REGEX REPLACE "," ";" _arch_list "${_target_architectures}")
    set(_all_accepted_architectures)
    foreach(_arch ${_arch_list})
        gmx_add_hipcc_flag_if_supported(GMX_HIP_HIPCC_FLAGS HIPCC_HAS_TARGET_ARCH_${_arch} "--offload-arch=${_arch}")
        if (HIPCC_HAS_TARGET_ARCH_${_arch})
            list(APPEND _all_accepted_architectures "${_arch}")
        endif()
    endforeach()
    if ("${_all_accepted_architectures}" STREQUAL "")
            message(FATAL_ERROR "No accepted offload target architectures were found for the hipcc compiler")
    endif()
    set(GMX_HIP_HIPCC_FLAGS ${GMX_HIP_HIPCC_FLAGS} PARENT_SCOPE)
    set(CMAKE_HIP_ARCHITECTURES ${_all_accepted_architectures} PARENT_SCOPE)
    # Only set the old AMDGPU_TARGETS for older versions of ROCm, use GPU_TARGETS for newer versions
    if (${HIP_VERSION} VERSION_LESS 6.4.0)
        set(AMDGPU_TARGETS ${_all_accepted_architectures} PARENT_SCOPE)
    else()
        set(GPU_TARGETS ${_all_accepted_architectures} PARENT_SCOPE)
    endif()
endfunction()

# iterate over user supplied and GROMACS default list of optimization flags
# to ensure that the compiler accepts them
function(gmx_hip_check_user_compile_flags _compile_flags)
    string(REGEX REPLACE "," ";" _compile_flags_list "${_compile_flags}")
    foreach(_flag ${_compile_flags_list})
        gmx_hip_check_single_flag(${_flag})
    endforeach()
    set(GMX_HIP_HIPCC_FLAGS ${GMX_HIP_HIPCC_FLAGS} PARENT_SCOPE)
endfunction()

# List of compilation flags used to control device and host code compilation. The functions used below and defined above all append to this list.
set(GMX_HIP_HIPCC_FLAGS)
# If the users doesn't limit the architectures to generate code for, we use a list of CDNA targets to build for by default.
# https://github.com/ROCm/ROCR-Runtime/blob/rocm-5.7.x/src/core/runtime/isa.cpp and https://rocm.docs.amd.com/_/downloads/HIP/en/latest/pdf/
set(GMX_HIP_TARGET_ARCH "gfx801,gfx802,gfx803,gfx900,gfx906,gfx90a,gfx90c,gfx942" CACHE STRING "Comma-separated list of target architectures to generate device code")

gmx_check_hip_architectures("${GMX_HIP_TARGET_ARCH}")

# Try to detect if we need RDNA support. Not very robust, but should cover the most common use.
if (${GMX_HIP_TARGET_ARCH} MATCHES "gfx1[0-9][0-9][0-9]")
    set(_enable_rdna_support_automatically ON)
else()
    set(_enable_rdna_support_automatically OFF)
    # We assume that any GCN2-5 architecture (gfx8) and CDNA1-3 (gfx9 series) up until the time of writing of this conditional is 64-wide
    if (${GMX_HIP_TARGET_ARCH} MATCHES "gfx8[0-9][0-9]|gfx9[0-4][0-9ac]")
        option(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
                "Disable NBNXM GPU cluster pair splitting. Only supported with HIP and 64-wide GPU architectures (like AMD GCN/CDNA)."
            ON)
        mark_as_advanced(GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT)
    endif()
endif()
option(GMX_ENABLE_AMD_RDNA_SUPPORT
        "Enable compiling kernels for AMD RDNA GPUs (gfx1xxx). When OFF, only CDNA and GCN are supported. Only used with HIP."
    ${_enable_rdna_support_automatically})
mark_as_advanced(GMX_ENABLE_AMD_RDNA_SUPPORT)

# Set the default compilation flags for building GROMACS with HIP. The defaults here can be appended to by changing either of the
# variables below, or by setting the GMX_HIPCC_EXTRA_FLAGS variable.
if (BUILD_SHARED_LIBS)
    gmx_hip_check_single_flag("-fPIC")
endif()
gmx_hip_check_single_flag("-fno-gpu-rdc")
gmx_hip_check_single_flag("-ffast-math")
gmx_hip_check_single_flag("-munsafe-fp-atomics")
gmx_hip_check_single_flag("-fdenormal-fp-math=ieee")
gmx_hip_check_single_flag("-fcuda-flush-denormals-to-zero")
gmx_hip_check_single_flag("-fno-slp-vectorize")
gmx_hip_check_single_flag("-Wno-unused-command-line-argument")
# currently ROCm 6.2.x spams warnings about missing occupancy targets during the
# compilation of the GPU kernels. We silence this warning here to prevent builds 
# failing because of that.
gmx_hip_check_single_flag("-Wno-pass-failed")

option(GMX_ENABLE_AMD_HIP_SAVE_TEMPS "Enable saving generated ISA to analyse compiler behaviour" FALSE)
mark_as_advanced(GMX_ENABLE_AMD_HIP_SAVE_TEMPS)
option(GMX_ENABLE_AMD_HIP_KERNEL_ANALYSIS "Enable generation of kernel analysis data" FALSE)
mark_as_advanced(GMX_ENABLE_AMD_HIP_KERNEL_ANALYSIS)

if (GMX_ENABLE_AMD_HIP_SAVE_TEMPS)
    gmx_hip_check_single_flag("--save-temps")
endif()

if (GMX_ENABLE_AMD_HIP_KERNEL_ANALYSIS)
    gmx_hip_check_single_flag("-Rpass-analysis=kernel-resource-usage")
endif()

# User may have supplied the optimization flags on the command line, only available for backwards compat with AMD port.
# In general we want to control those flags based on the GROMACS build type, but won't stop users supplying different
# options if they think they know better.
gmx_hip_check_user_compile_flags("${GMX_HIPCC_EXTRA_FLAGS}")

# In case any additional flags have been set to control the HIPCC compiler, we add them to our list of compilation flags.
gmx_hip_check_user_compile_flags("${HIPCC_EXTRA_FLAGS}")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
if(${_cmake_build_type} STREQUAL "RELEASE")
    gmx_hip_check_single_flag("-DNDEBUG")
elseif(${_cmake_build_type} MATCHES "DEB") # should match all builds that include debug symbols, but exclude other RELWITH* builds
    gmx_hip_check_single_flag("-ggdb")
    gmx_hip_check_single_flag("-O0")
endif()

# This helper function creates a temporary scope in which we can set
# the definitions, include directories and add the correct compiler flags
# and properties to allow correct compilation with the HIP compiler.
# We are not using the wrapper that comes with the ROCm toolkit but are
# relying directly on the built-in support in CMake to do that.
function(gmx_hip_add_library TARGET)
    add_definitions(-DHAVE_CONFIG_H)

    # Now add all the compilation options
    gmx_device_target_compile_options(HIP_${TARGET}_CXXFLAGS)
    gmx_hip_check_user_compile_flags("${HIP_${TARGET}_CXXFLAGS}")
    gmx_hip_check_user_compile_flags("${HIP_${TARGET}_CXXFLAGS_${CMAKE_BUILD_TYPE}}")

    add_library(${TARGET} ${ARGN})
    set_property(TARGET ${TARGET} PROPERTY HIP_STANDARD ${CMAKE_CXX_STANDARD})
    set_target_properties(${TARGET} PROPERTIES HIP_ARCHITECTURES OFF)
    target_compile_options(${TARGET} PRIVATE $<$<COMPILE_LANGUAGE:HIP>:${GMX_HIP_HIPCC_FLAGS}>)

    # TODO: Restrict the scope of MPI dependence.
    # Targets that actually need MPI headers and build tool flags should
    # manage their own `target_link_libraries` locally. Such a change is beyond
    # the scope of the bug fix for #4678.
    if (GMX_LIB_MPI)
        target_link_libraries(${TARGET} PRIVATE MPI::MPI_CXX)
    endif ()
endfunction()
