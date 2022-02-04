#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
# Copyright (c) 2017,2018,2019,2020,2021,2022, by the GROMACS development team, led by
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

set(GMX_GPU_CUDA ON)

option(GMX_CLANG_CUDA "Use clang for CUDA" OFF)

if(GMX_DOUBLE)
    message(FATAL_ERROR "CUDA acceleration is not available in double precision")
endif()

find_package(CUDA ${REQUIRED_CUDA_VERSION} REQUIRED)

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

option(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT "Whether to compile the CUDA non-bonded module using a single compilation unit." OFF)
mark_as_advanced(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT)
