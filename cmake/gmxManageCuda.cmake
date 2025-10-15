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

set(GMX_GPU_CUDA ON)

option(GMX_CLANG_CUDA "Use clang for CUDA" OFF)

if(GMX_DOUBLE)
    message(FATAL_ERROR "CUDA acceleration is not available in double precision")
endif()

set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

find_package(CUDAToolkit ${GMX_CUDA_MINIMUM_REQUIRED_VERSION} REQUIRED)
set(GMX_HAVE_GPU_GRAPH_SUPPORT ON)

if(NOT CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
endif()
if(GMX_CLANG_CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        set(CMAKE_CUDA_COMPILER "${CMAKE_CXX_COMPILER}")
    endif()
else()
    # Using NVIDIA compiler
    if(NOT CUDAToolkit_NVCC_EXECUTABLE)
        message(FATAL_ERROR "nvcc is required for a CUDA build, please set CUDAToolkit_ROOT appropriately")
    endif()
    if (NOT CMAKE_CUDA_COMPILER)
        set(CMAKE_CUDA_COMPILER "${CUDAToolkit_NVCC_EXECUTABLE}")
    endif()
endif()

# We do this before enable_language(CUDA) to avoid CMAKE_CUDA_ARCHITECTURES being set to whatever enable_language(CUDA) autodetects
# But if we _do_ want to use the autodetected architectures, we move enable_language(CUDA) before this block
if(NOT GMX_CUDA_ARCHITECTURES)
    if (GMX_CUDA_TARGET_SM OR GMX_CUDA_TARGET_COMPUTE)
        message(WARNING "GMX_CUDA_TARGET_SM and GMX_CUDA_TARGET_COMPUTE are deprecated, use CMAKE_CUDA_ARCHITECTURES instead")
        set(GMX_CUDA_ARCHITECTURES ${GMX_CUDA_TARGET_SM})
        set(_target_compute_list ${GMX_CUDA_TARGET_COMPUTE})
        foreach(_target ${_target_compute_list})
            list(APPEND GMX_CUDA_ARCHITECTURES "${_target}-virtual")
        endforeach()
    elseif(NOT CMAKE_CUDA_ARCHITECTURES)
        set(GMX_CUDA_ARCHITECTURES "all")
    else()
        set(GMX_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
    endif()
endif()

enable_language(CUDA)

# CMAKE_CUDA_ARCHITECTURES_ALL is set by enable_language(CUDA)
# We expand the list so that it can be pruned down later if NVSHMEM is enabled
if(GMX_CUDA_ARCHITECTURES STREQUAL "all")
    set(GMX_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES_ALL}")
endif()

set(GMX_CUDA_ARCHITECTURES "${GMX_CUDA_ARCHITECTURES}" CACHE STRING "GROMACS-specific CUDA architectures")
mark_as_advanced(GMX_CUDA_ARCHITECTURES)

set(_cuda_arch_message "${GMX_CUDA_ARCHITECTURES}")
string(SHA1 _cuda_arch_message_hash "${_cuda_arch_message}")
if(NOT LAST_REPORTED_CUDA_ARCH_MESSAGE_HASH_CUDA STREQUAL _cuda_arch_message_hash)
    message(STATUS "Compiling GROMACS for CUDA architectures: ${_cuda_arch_message}")
    set(LAST_REPORTED_CUDA_ARCH_MESSAGE_HASH_CUDA "${_cuda_arch_message_hash}" CACHE INTERNAL "The hash of the last reported CUDA architecture list, after the initial CUDA detection")
endif()

if (GMX_OPENMP AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.31)
    # Re-run OpenMP detection to add OpenMP::OpenMP_CUDA target
    if (OpenMP_CUDA_FIND_QUIETLY_AFTER_FIRST_RUN AND OpenMP_CXX_FIND_QUIETLY_AFTER_FIRST_RUN)
        set (OpenMP_FIND_QUIETLY TRUE)
    else()
        set (OpenMP_FIND_QUIETLY FALSE)
    endif()
    find_package(OpenMP REQUIRED COMPONENTS CXX CUDA)
    set(OpenMP_CXX_FIND_QUIETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to find OpenMP_CXX")
    set(OpenMP_CUDA_FIND_QUIETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to find OpenMP_CUDA")
endif()

# Tests a single flag to use with CUDA compiler.
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
include(CheckCompilerFlag)
function(gmx_add_cuda_flag_if_supported _flags_cache_variable_name)
    check_compiler_flag(CUDA "${ARGN}" ${_flags_cache_variable_name})
    if (${_flags_cache_variable_name})
        list(APPEND GMX_CUDA_FLAGS ${ARGN})
        set(GMX_CUDA_FLAGS "${GMX_CUDA_FLAGS}" PARENT_SCOPE)
    endif()
endfunction()

set(GMX_CUDA_FLAGS)
if(GMX_CLANG_CUDA)
    include(gmxManageClangCudaConfig)
else()
    include(gmxManageNvccConfig)
endif()
set(GMX_CUDA_FLAGS "${GMX_CUDA_FLAGS}" CACHE STRING "GROMACS-specific CUDA flags")
mark_as_advanced(GMX_CUDA_FLAGS)

option(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT "Whether to compile the CUDA non-bonded module using a single compilation unit." OFF)
mark_as_advanced(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT)

macro(get_cuda_compiler_info COMPILER_INFO COMPILER_FLAGS COMPILER_ARCHITECTURES)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    # CXX compiler is the CUDA compiler
    set(${COMPILER_INFO} "${CMAKE_CUDA_COMPILER} (${CMAKE_CUDA_COMPILER_ID} ${CMAKE_CUDA_COMPILER_VERSION})")
    # there are some extra flags
    list(JOIN GMX_CUDA_FLAGS " " GMX_CUDA_FLAGS_STR)
    set(${COMPILER_FLAGS} "${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${_cmake_build_type}} ${GMX_CUDA_FLAGS_STR}")
    set(${COMPILER_ARCHITECTURES} ${GMX_CUDA_ARCHITECTURES})
endmacro ()
