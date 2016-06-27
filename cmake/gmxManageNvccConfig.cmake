#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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
# - set icc compatibility mode to gcc 4.6
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

# Set the host compiler for nvcc if this is not set by CMake (v<=2.8.9)
#
# Note that even though nvcc compiles host code as C++, we use the
# CMAKE_C_COMPILER as host compiler. We do this because CUDA versions
# preceding 5.0 only recognize icc, but not icpc. However, both gcc and icc
# (i.e. all supported compilers) happily compile C++ code.
#
# Also note that with MSVC nvcc sets the -compiler-bindir option behind the
# scenes; to avoid conflicts we don't set -ccbin automatically.
#
# TODO: remove this when CMAke >=v2.8.10 is required.
if (NOT DEFINED CUDA_HOST_COMPILER AND NOT MSVC)
    set(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}")
    set(CUDA_HOST_COMPILER_AUTOSET TRUE CACHE INTERNAL
        "True if CUDA_HOST_COMPILER is automatically set")
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
    # FindCUDA in CMake 2.8.10 sets the host compiler internally
    if (CMAKE_VERSION VERSION_LESS "2.8.10")
        set(CUDA_HOST_COMPILER ${CUDA_HOST_COMPILER}
            CACHE PATH "Host compiler for nvcc")
    endif()

    # On *nix force icc in gcc 4.6 compatibility mode. This is needed
    # as even with icc used as host compiler, when icc's gcc compatibility
    # mode is higher than the max gcc version officially supported by CUDA,
    # nvcc will freak out.
    set(CUDA_HOST_COMPILER_OPTIONS "")
    if (UNIX AND
            ((CMAKE_C_COMPILER_ID MATCHES "Intel" AND
              (CUDA_HOST_COMPILER_AUTOSET OR CMAKE_C_COMPILER STREQUAL CUDA_HOST_COMPILER)) OR
            (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER STREQUAL CUDA_HOST_COMPILER))
        )
        message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.6 for nvcc host compilation")
        list(APPEND CUDA_HOST_COMPILER_OPTIONS "-Xcompiler;-gcc-version=460")
    endif()

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
    #
    #
    #   Note that CUDA 6.5.19 second patch release supports cc 5.2 too, but
    #   CUDA_VERSION does not contain patch version and having PTX 5.0 JIT-ed is
    #   equally fast as compiling with sm_5.2 anyway.

    # First add flags that trigger SASS (binary) code generation for physical arch
    list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_20,code=sm_20")
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

    # Next add flags that trigger PTX code generation for the newest supported virtual arch
    # that's useful to JIT to future architectures
    if(CUDA_VERSION VERSION_LESS "6.5")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_35,code=compute_35")
    elseif(CUDA_VERSION VERSION_LESS "7.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_50,code=compute_50")
    elseif(CUDA_VERSION VERSION_LESS "8.0")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_52,code=compute_52")
    else() # version >= 8.0
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_60,code=compute_60")
        list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_61,code=compute_61")
    endif()
endif()

gmx_dependent_cache_variable(GMX_CUDA_TARGET_SM "List of CUDA GPU architecture codes to compile for (without the sm_ prefix)" STRING "" GMX_CUDA_TARGET_SM)
gmx_dependent_cache_variable(GMX_CUDA_TARGET_COMPUTE "List of CUDA virtual architecture codes to compile for (without the compute_ prefix)" STRING "" GMX_CUDA_TARGET_COMPUTE)

# assemble the CUDA flags
list(APPEND GMX_CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_GENCODE_FLAGS}")
list(APPEND GMX_CUDA_NVCC_FLAGS "-use_fast_math")

# assemble the CUDA host compiler flags
# with CMake <2.8.10 the host compiler needs to be set on the nvcc command line
if (CMAKE_VERSION VERSION_LESS "2.8.10")
    list(APPEND GMX_CUDA_NVCC_FLAGS "-ccbin=${CUDA_HOST_COMPILER}")
endif()
list(APPEND GMX_CUDA_NVCC_FLAGS "${CUDA_HOST_COMPILER_OPTIONS}")

# The flags are set as local variables which shadow the cache variables. The cache variables
# (can be set by the user) are appended. This is done in a macro to set the flags when all
# host compiler flags are already set.
macro(GMX_SET_CUDA_NVCC_FLAGS)
    if(CUDA_PROPAGATE_HOST_FLAGS)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)

        # When CUDA 6.5 is required we should use C++11 also for CUDA and also propagate
        # the C++11 flag to CUDA. Then we can use the solution implemented in FindCUDA
        # (starting with 3.3 - can be backported). For now we need to remove the C++11
        # flag which means we need to manually propagate all other flags.
        string(REGEX REPLACE "[-]+std=c\\+\\+0x" "" _CMAKE_CXX_FLAGS_SANITIZED "${CMAKE_CXX_FLAGS}")

        # The IBM xlc compiler chokes if we use both altivec and Cuda. Solve
        # this by not propagating the flag in this case.
        if(CMAKE_CXX_COMPILER_ID MATCHES "XL")
            string(REGEX REPLACE "-qaltivec" "" _CMAKE_CXX_FLAGS_SANITIZED "${_CMAKE_CXX_FLAGS_SANITIZED}")
        endif()

        # CUDA versions prior to 7.5 come with a header (math_functions.h) which uses the _MSC_VER macro
        # unconditionally, so we strip -Wundef from the propagatest flags for earlier CUDA versions.
        if (CUDA_VERSION VERSION_LESS "7.5")
            string(REGEX REPLACE "-Wundef" "" _CMAKE_CXX_FLAGS_SANITIZED "${_CMAKE_CXX_FLAGS_SANITIZED}")
        endif()

        string(REPLACE " " "," _flags "${_CMAKE_CXX_FLAGS_SANITIZED}")
        set(CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_FLAGS};${CUDA_NVCC_FLAGS};-Xcompiler;${_flags}")

        # Create list of all possible configurations. For multi-configuration this is CMAKE_CONFIGURATION_TYPES
        # and for single configuration CMAKE_BUILD_TYPE. Not sure why to add the default ones, but FindCUDA
        # claims one should.
        set(CUDA_configuration_types ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE} Debug MinSizeRel Release RelWithDebInfo)
        list(REMOVE_DUPLICATES CUDA_configuration_types)

        foreach(_config ${CUDA_configuration_types})
            string(TOUPPER ${_config} _config_upper)
            string(REPLACE " " "," _flags "${CMAKE_CXX_FLAGS_${_config_upper}}")
            set(CUDA_NVCC_FLAGS_${_config_upper} "${CUDA_NVCC_FLAGS_${_config_upper}};-Xcompiler;${_flags}")
        endforeach()
    else()
        set(CUDA_NVCC_FLAGS "${GMX_CUDA_NVCC_FLAGS};${CUDA_NVCC_FLAGS}")
    endif()
endmacro()
