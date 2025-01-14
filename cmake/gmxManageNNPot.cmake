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
    "Enable neural network potential interface."
    AUTO
    AUTO TORCH OFF
)

if(TORCH_ALREADY_SEARCHED)
    set(FIND_TORCH_QUIETLY ON)
endif()

if(NOT GMX_NNPOT STREQUAL "OFF")

    find_package(Torch 2.0.0 QUIET)
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