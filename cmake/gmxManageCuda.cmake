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

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

find_package(CUDA ${REQUIRED_CUDA_VERSION} REQUIRED)

if(${CUDA_VERSION} GREATER_EQUAL 11.1)
  set(GMX_HAVE_GPU_GRAPH_SUPPORT ON)
endif()

mark_as_advanced(CUDA_SDK_ROOT_DIR CUDA_USE_STATIC_CUDA_RUNTIME)

# Try to execute ${CUDA_NVCC_EXECUTABLE} --version and set the output
# (or an error string) in the argument variable.
# Note that semicolon is used as separator for nvcc.
#
# Parameters:
#   COMPILER_INFO         - [output variable] string with compiler path, ID and
#                           some compiler-provided information
#   DEVICE_COMPILER_FLAGS - [output variable] device flags for the compiler
#   HOST_COMPILER_FLAGS   - [output variable] host flags for the compiler, if propagated
#
macro(get_cuda_compiler_info COMPILER_INFO DEVICE_COMPILER_FLAGS HOST_COMPILER_FLAGS)
    if(NOT GMX_CLANG_CUDA)
        if(CUDA_NVCC_EXECUTABLE)

            # Get the nvcc version string. This is multi-line, but since it is only 4 lines
            # and might change in the future it is better to store than trying to parse out
            # the version from the current format.
            execute_process(COMMAND ${CUDA_NVCC_EXECUTABLE} --version
                RESULT_VARIABLE _nvcc_version_res
                OUTPUT_VARIABLE _nvcc_version_out
                ERROR_VARIABLE  _nvcc_version_err
                OUTPUT_STRIP_TRAILING_WHITESPACE)
            if (${_nvcc_version_res} EQUAL 0)
                # Fix multi-line mess: Replace newline with ";" so we can use it in a define
                string(REPLACE "\n" ";" _nvcc_info_singleline ${_nvcc_version_out})
                SET(${COMPILER_INFO} "${CUDA_NVCC_EXECUTABLE} ${_nvcc_info_singleline}")
                string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
                SET(_compiler_flags "${CUDA_NVCC_FLAGS_${_build_type}}")
                if(CUDA_PROPAGATE_HOST_FLAGS)
                    set(${HOST_COMPILER_FLAGS} BUILD_CXXFLAGS)
                else()
                    set(${HOST_COMPILER_FLAGS} "")
                endif()
                SET(${DEVICE_COMPILER_FLAGS} "${CUDA_NVCC_FLAGS}${CUDA_NVCC_FLAGS_${_build_type}}")
            else()
                SET(${COMPILER_INFO} "N/A")
                SET(${COMPILER_FLAGS} "N/A")
            endif()
        endif()
    else()
        # CXX compiler is the CUDA compiler
        set(${COMPILER_INFO} "${CMAKE_CXX_COMPILER}  ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
        # there are some extra flags
        set(${COMPILER_FLAGS} "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}} ${GMX_CUDA_CLANG_FLAGS}")
    endif()
endmacro ()

if(GMX_CLANG_CUDA)
    include(gmxManageClangCudaConfig)
    list(APPEND GMX_EXTRA_LIBRARIES ${GMX_CUDA_CLANG_LINK_LIBS})
    link_directories("${GMX_CUDA_CLANG_LINK_DIRS}")
else()
    # Using NVIDIA compiler
    if(NOT CUDA_NVCC_EXECUTABLE)
        message(FATAL_ERROR "nvcc is required for a CUDA build, please set CUDA_TOOLKIT_ROOT_DIR appropriately")
    endif()
    # set up nvcc options
    include(gmxManageNvccConfig)
endif()

# We make sure to call get_cuda_compiler_info() before reaching this line,
# so we report errors related to host compiler / nvcc mismatch
# before the call to enable_language(CUDA).
enable_language(CUDA)

option(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT "Whether to compile the CUDA non-bonded module using a single compilation unit." OFF)
mark_as_advanced(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT)

# custom libcuda and libnvdia-ml stub lib finder as find_package(CUDA) doesn't support it.
if (GMX_USE_CUFFTMP OR GMX_NVSHMEM)
    find_library(GMX_CUDA_DRV_LIB cuda HINTS "${CUDA_TOOLKIT_ROOT_DIR}" PATH_SUFFIXES "lib64/stubs" REQUIRED)
    find_library(GMX_NVIDIA_ML_LIB nvidia-ml HINTS "${CUDA_TOOLKIT_ROOT_DIR}" PATH_SUFFIXES "lib64/stubs" REQUIRED)
endif()
