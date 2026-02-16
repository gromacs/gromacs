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

# Add a flag to enable neural network potential support. Currently only supports libtorch.
gmx_option_multichoice(GMX_NNPOT
    "Enable neural network potential interface"
    AUTO
    AUTO TORCH OFF
)

if(TORCH_ALREADY_SEARCHED)
    set(FIND_TORCH_QUIETLY ON)
endif()

set(GMX_TORCH OFF)
if(NOT GMX_NNPOT STREQUAL "OFF")
    if(GMX_GPU_CUDA AND NOT TORCH_CUDA_ARCH_LIST)
        set(TORCH_CUDA_ARCH_LIST)
        foreach(_arch IN LISTS GMX_CUDA_ARCHITECTURES)
            if(_arch MATCHES "^[0-9]+[a-z]?(-virtual)?$")
                # Convert _arch from 75 or 75-virtual to 7.5+PTX
                string(REGEX REPLACE "^([0-9]+)([0-9][a-z]?)(-virtual)?$" "\\1.\\2+PTX" arch_ptx "${_arch}")
            elseif(_arch MATCHES "^[0-9]+[a-z]?-real$")
                # Convert _arch from 75-real to 7.5
                string(REGEX REPLACE "^([0-9]+)([0-9][a-z]?)-real$" "\\1.\\2" arch_ptx "${_arch}")
            else()
                message(FATAL_ERROR "Unknown CUDA architecture: ${_arch}")
            endif()
            set(TORCH_CUDA_ARCH_LIST "${TORCH_CUDA_ARCH_LIST} ${arch_ptx}")
        endforeach()
    endif()

    # Not so nice workaround because torch disables CMAKE_CUDA_ARCHITECTURES and complains if it is set
    set(_cmake_cuda_architectures_bak "${CMAKE_CUDA_ARCHITECTURES}")
    unset(CMAKE_CUDA_ARCHITECTURES CACHE) # yikes!
    # When we require at least CMake 4.1, finding Torch with OPTIONAL might be a good approach
    find_package(Torch 2.0.0 QUIET)
    set(CMAKE_CUDA_ARCHITECTURES "${_cmake_cuda_architectures_bak}" CACHE STRING "")
    set(TORCH_ALREADY_SEARCHED TRUE CACHE BOOL "True if a search for libtorch has already been done")
    mark_as_advanced(TORCH_ALREADY_SEARCHED)

    if(Torch_FOUND)
        # TORCH_LIBRARIES contain imported target "torch" that will set all flags and include paths etc
        list(APPEND GMX_COMMON_LIBRARIES ${TORCH_LIBRARIES})
        if(NOT FIND_TORCH_QUIETLY)
            message(STATUS "Found Torch: Neural network potential support enabled.")
        endif()

        # Check if the Torch version uses the correct ABI
        if (${TORCH_CXX_FLAGS} MATCHES "-D_GLIBCXX_USE_CXX11_ABI=0")
            message(FATAL_ERROR "Torch was compiled with the pre-cxx11 ABI. Please use a libtorch version "
                                "compiled with the cxx11 ABI, which is required for building GROMACS.")
        endif()

        set(GMX_TORCH ON)
    elseif(GMX_NNPOT STREQUAL "TORCH")
        message(FATAL_ERROR "Torch not found. Please install libtorch and add its installation prefix"
                            " to CMAKE_PREFIX_PATH or set Torch_DIR to a directory containing "
                            "a TorchConfig.cmake or torch-config.cmake file.")
    else() # "AUTO"
        if(NOT FIND_TORCH_QUIETLY)
            message(STATUS "Torch not found. Neural network potential support will be disabled.")
        endif()
    endif()
endif()
