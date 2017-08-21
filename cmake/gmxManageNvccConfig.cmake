#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
# - auto-detect compatible nvcc host compiler and set nvcc -ccbin (if not MPI wrapper)
# - set icc compatibility mode to gcc 4.8.1
# - (advanced) variables set:
#   * CUDA_HOST_COMPILER            - the host compiler for nvcc (only with cmake <2.8.10)
#   * CUDA_HOST_COMPILER_OPTIONS    - the full host-compiler related option list passed to nvcc
#
# Note that from CMake 2.8.10 FindCUDA defines CUDA_HOST_COMPILER internally,
# so we won't set it ourselves, but hope that the module does a good job.

gmx_check_if_changed(CUDA_HOST_COMPILER_CHANGED CUDA_HOST_COMPILER)

# CUDA_HOST_COMPILER changed hence it is not auto-set anymore
if (CUDA_HOST_COMPILER_CHANGED AND CUDA_HOST_COMPILER_AUTOSET)
    unset(CUDA_HOST_COMPILER_AUTOSET CACHE)
endif()

# glibc 2.23 changed string.h in a way that breaks CUDA compilation in
# many projects, but which has a trivial workaround. It would be nicer
# to compile with nvcc and see that the workaround is necessary and
# effective, but it is unclear how to do that. Also, grepping in the
# glibc source shows that _FORCE_INLINES is only used in this string.h
# feature and performance of memcpy variants is unimportant for CUDA
# code in GROMACS. So this workaround is good enough to keep problems
# away from users installing GROMACS. See Redmine 1942.
function(work_around_glibc_2_23)
    try_compile(IS_GLIBC_2_23_OR_HIGHER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/TestGlibcVersion.cpp)
    if(IS_GLIBC_2_23_OR_HIGHER)
        message(STATUS "Adding work-around for issue compiling CUDA code with glibc 2.23 string.h")
        list(APPEND CUDA_HOST_COMPILER_OPTIONS "-D_FORCE_INLINES")
        set(CUDA_HOST_COMPILER_OPTIONS ${CUDA_HOST_COMPILER_OPTIONS} PARENT_SCOPE)
    endif()
endfunction()

# set up host compiler and its options
if(CUDA_HOST_COMPILER_CHANGED)
    set(CUDA_HOST_COMPILER_OPTIONS "")

    if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "GNU")
        # Some versions of gcc-4.8 and gcc-4.9 produce errors (in particular on OS X)
        # if we do not use -D__STRICT_ANSI__. It is harmless, so we might as well add it for all versions.
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
    # - with CUDA >=5.0 <6.5:   CC <=3.5 is supported
    #     => compile sm_20, sm_30, sm_35 SASS, and compute_35 PTX
    # - with CUDA ==6.5:        CC <=3.7 and 5.0 are supported
    #     => compile sm_20, sm_30, sm_35, sm_37 sm_50, SASS, and compute_50 PTX
    # - with CUDA >=7.0         CC 5.2 is supported (5.3, Tegra X1 we don't generate code for)
    #     => compile sm_20, sm_30, sm_35, sm_37, sm_50, & sm_52 SASS, and compute_52 PTX
    # - with CUDA >=8.0         CC 6.0-6.2 is supported (but we know nothing about CC 6.2, so we won't generate code or it)
    #     => compile sm_20, sm_30, sm_35, sm_37, sm_50, sm_52, sm_60, sm_61 SASS, and compute_60 and compute_61 PTX
    # - with CUDA >=9.0         CC 7.0 is supported and CC 2.0 is no longer supported
    #     => compile sm_30, sm_35, sm_37, sm_50, sm_52, sm_60, sm_61, sm_70 SASS, and compute_70 PTX
    #
    #   Note that CUDA 6.5.19 second patch release supports cc 5.2 too, but
    #   CUDA_VERSION does not contain patch version and having PTX 5.0 JIT-ed is
    #   equally fast as compiling with sm_5.2 anyway.

    # First add flags that trigger SASS (binary) code generation for physical arch
    if(CUDA_VERSION VERSION_LESS "9.00") # < 9.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_20,code=sm_20")
    endif()
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_30,code=sm_30")
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=sm_35")

    if(NOT CUDA_VERSION VERSION_LESS "6.5") # >= 6.5
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_37,code=sm_37")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=sm_50")
    endif()
    if(NOT CUDA_VERSION VERSION_LESS "7.0") # >= 7.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_52,code=sm_52")
    endif()
    if(NOT CUDA_VERSION VERSION_LESS "8.0") # >= 8.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_60,code=sm_60")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_61,code=sm_61")
    endif()
    if(NOT CUDA_VERSION VERSION_LESS "9.0") # >= 9.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_70,code=sm_70")
    endif()

    # Next add flags that trigger PTX code generation for the newest supported virtual arch
    # that's useful to JIT to future architectures
    if(CUDA_VERSION VERSION_LESS "6.5")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=compute_35")
    elseif(CUDA_VERSION VERSION_LESS "7.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=compute_50")
    elseif(CUDA_VERSION VERSION_LESS "8.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_52,code=compute_52")
    elseif(CUDA_VERSION VERSION_LESS "9.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_60,code=compute_60")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_61,code=compute_61")
    else() # version >= 9.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_70,code=compute_70")
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
if (CUDA_VERSION VERSION_EQUAL "8.0")
    # requesting sm_20 triggers deprecation messages with nvcc 8.0 which we better avoid
    list(APPEND GMX_CUDA_NVCC_FLAGS "-Wno-deprecated-gpu-targets")
endif()

# assemble the CUDA host compiler flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${CUDA_HOST_COMPILER_OPTIONS}")

# The flags are set as local variables which shadow the cache variables. The cache variables
# (can be set by the user) are appended. This is done in a macro to set the flags when all
# host compiler flags are already set.
macro(GMX_SET_CUDA_NVCC_FLAGS)
    set(CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_FLAGS};${CUDA_NVCC_FLAGS}")
endmacro()
