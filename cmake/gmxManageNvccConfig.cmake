#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,215, by the GROMACS development team, led by
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
# - set icc compatibility mode to gcc 4.4/4.5 (CUDA 4.0 is not compatible with gcc >v4.4)
# - (advanced) variables set:
#   * CUDA_HOST_COMPILER_OPTIONS    - the full host-compiler related option list passed to nvcc
#

gmx_check_if_changed(CUDA_HOST_COMPILER_CHANGED CUDA_HOST_COMPILER)

# set up host compiler and its options
if(CUDA_HOST_COMPILER_CHANGED)
    # On *nix force icc in gcc 4.4 compatibility mode with CUDA 3.2/4.0 and
    # gcc 4.5 compatibility mode with later CUDA versions. This is needed
    # as even with icc used as host compiler, when icc's gcc compatibility
    # mode is higher than the max gcc version officially supported by CUDA,
    # nvcc will freak out.
    set(CUDA_HOST_COMPILER_OPTIONS "")
    if (UNIX AND
            ((CMAKE_C_COMPILER_ID MATCHES "Intel" AND
              CMAKE_C_COMPILER STREQUAL CUDA_HOST_COMPILER) OR
            (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER STREQUAL CUDA_HOST_COMPILER))
        )
        if (CUDA_VERSION VERSION_LESS "4.1")
            message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.4 for nvcc host compilation")
            list(APPEND CUDA_HOST_COMPILER_OPTIONS "-Xcompiler;-gcc-version=440")
        else()
            message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.5 for nvcc host compilation")
            list(APPEND CUDA_HOST_COMPILER_OPTIONS "-Xcompiler;-gcc-version=450")
        endif()
    endif()

    if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "GNU")
        # Some versions of gcc-4.8 and gcc-4.9 produce errors (in particular on OS X)
        # if we do not use -D__STRICT_ANSI__. It is harmless, so we might as well add it for all versions.
        list(APPEND CUDA_HOST_COMPILER_OPTIONS "-D__STRICT_ANSI__")
    endif()

    set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS}"
        CACHE STRING "Options for nvcc host compiler (do not edit!).")

    mark_as_advanced(CUDA_HOST_COMPILER_OPTIONS)
endif()

# the legacy CUDA kernels have been dropped, warn with CUDA 4.0
if (CUDA_VERSION VERSION_EQUAL "4.0")
    message(WARNING "The legacy GPU kernels optimized for older CUDA compilers, including the detected version 4.0, have been removed. To avoid performance loss, we strongly recommend upgrading to a newer CUDA toolkit.
    ")
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
    # - with CUDA  <4.2:        compute capability 2.x supported (compiling for sm_2.1 does not help):
    #     => compile sm_20 cubin, and compute_20 PTX
    # - with CUDA  =4.2 <5.0:   CC <=3.0 is supported:
    #     => compile sm_20, sm_30 cubin, and compute_30 PTX
    # - with CUDA >=5.0 <6.5:   CC <=3.5 is supported
    #     => compile sm_20, sm_30, sm_35 cubin, and compute_35 PTX
    # - with CUDA >=6.5:        CC <=3.7 and 5.0 are supported
    #     => compile sm_20, sm_30, sm_35, sm_37 sm_50, cubin, and compute_50 PTX
    #
    #   Note that CUDA 6.5.19 second patch release supports cc 5.2 too, but
    #   CUDA_VERSION does not contain patch version and having PTX 5.0 JIT-ed is
    #   equally fast as compiling with sm_5.2 anyway.
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_20,code=sm_20")
    if(CUDA_VERSION VERSION_GREATER "4.1990") # >= 4.2
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_30,code=sm_30")
    endif()
    if(CUDA_VERSION VERSION_GREATER "4.999") # >= 5.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=sm_35")
    endif()
    if(CUDA_VERSION VERSION_GREATER "6.4999") # >= 6.5
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_37,code=sm_37")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=sm_50")
    endif()

    if(CUDA_VERSION VERSION_LESS "4.2")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_20,code=compute_20")
    elseif(CUDA_VERSION VERSION_LESS "5.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_30,code=compute_30")
    elseif(CUDA_VERSION VERSION_LESS "6.5")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=compute_35")
    else() # version >= 6.5
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=compute_50")
    endif()
endif()

gmx_dependent_cache_variable(GMX_CUDA_TARGET_SM "List of CUDA GPU architecture codes to compile for (without the sm_ prefix)" STRING "" GMX_CUDA_TARGET_SM)
gmx_dependent_cache_variable(GMX_CUDA_TARGET_COMPUTE "List of CUDA virtual architecture codes to compile for (without the compute_ prefix)" STRING "" GMX_CUDA_TARGET_COMPUTE)

# assemble the CUDA flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_GENCODE_FLAGS}")
list(APPEND GMX_CUDA_NVCC_FLAGS "-use_fast_math")

# assemble the CUDA host compiler flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${CUDA_HOST_COMPILER_OPTIONS}")

# finally set the damn flags
set(CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_FLAGS}" CACHE STRING "Compiler flags for nvcc." FORCE)
