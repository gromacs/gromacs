#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
# - set icc compatibility mode to gcc 4.4/4.5 (CUDA 4.0 is not compatible with gcc >v4.4)
# - (advanced) variables set:
#   * CUDA_HOST_COMPILER            - the host compiler for nvcc (only with cmake <2.8.10)
#   * CUDA_HOST_COMPILER_OPTIONS    - the full host-compiler related option list passed to nvcc
#
# Note that from CMake 2.8.10 FindCUDA defines CUDA_HOST_COMPILER internally,
# so we won't set it ourselves, but hope that the module does a good job.

if (NOT DEFINED CUDA_NVCC_FLAGS_SET)
    set(CUDA_NVCC_FLAGS_SET TRUE CACHE INTERNAL "True if NVCC flags have been set" FORCE)

    # Explicitly set the host compiler for nvcc if the current compiler is
    # supported and it's not an MPI compiler wrapper, otherwise warn the user.
    #
    # Note that even though nvcc compiles host code as C++, we use the
    # CMAKE_C_COMPILER as host compiler. We do this because CUDA versions
    # preceding 5.0 only recognize icc, but not icpc. However, both gcc and icc
    # (i.e. all supported compilers) happily compile C++ code.
    #
    # Also note that with MSVC nvcc sets the -compiler-bindir option behind the
    # scenes; to avoid conflicts we don't set -ccbin automatically.

    if (NOT DEFINED CUDA_HOST_COMPILER AND NOT MSVC)
        if (NOT CMAKE_COMPILER_IS_GNUCC AND
            NOT (CMAKE_C_COMPILER_ID MATCHES "Intel" AND UNIX AND NOT APPLE))
            message(WARNING "
            Will not set the nvcc host compiler because the current C compiler is not
            compatible with nvcc:
            ${CMAKE_C_COMPILER} (ID: ${CMAKE_C_COMPILER_ID})
            Compatible compilers are: gcc on Linux and Mac OS X, the Intel Compiler on 64-bit
            Linux and MSVC on Windows. Note that with newer CUDA releases this might change,
            for up-to-date compatibility information check the NVIDIA documentation.
            If nothing specified, nvcc will automatically pick the platform-default compiler;
            Note that mixing compilers can cause errors.
            To manually set the nvcc host compiler, edit CUDA_NVCC_FLAGS or re-configure
            setting CUDA_HOST_COMPILER to the full path of a compatible compiler.
            ")
        else()
            # do not use MPI compiler wrappers, as these are prone to brake nvcc
            if (GMX_MPI AND
                NOT "${${MPI_PREFIX}_FOUND}" AND # FindMPI-based detection
                NOT GMX_THREAD_MPI)
                message(WARNING "
            Will not set the nvcc host compiler because the current C compiler is an MPI
            compiler wrapper: ${CMAKE_C_COMPILER}
            MPI compiler wrappers are prone to not work with nvcc. You might get lucky,
            but the safest is to use the C compiler that the MPI compiler wrapper uses
            (if this is compatible).
            To manually set the nvcc host compiler, edit CUDA_NVCC_FLAGS or re-configure
            setting CUDA_HOST_COMPILER to the full path of a compatible compiler.
            ")
            else()
                set(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}")
                set(CUDA_HOST_COMPILER_AUTOSET TRUE CACHE INTERNAL
                    "True if CUDA_HOST_COMPILER is automatically set" FORCE)
            endif()
        endif()
    endif()

    if(DEFINED CUDA_HOST_COMPILER)
        # FindCUDA in CMake 2.8.10 sets the host compiler internally
        if (CMAKE_VERSION VERSION_LESS "2.8.10")
            message(STATUS "Setting the nvcc host compiler to: ${CUDA_HOST_COMPILER}")
            set(CUDA_HOST_COMPILER ${CUDA_HOST_COMPILER}
                CACHE PATH "Host compiler for nvcc (do not edit!)" FORCE)
            set(_HOST_COMPILER_OPTION_STRING "-ccbin=${CUDA_HOST_COMPILER};")
        endif()

        # On *nix force icc in gcc 4.4 compatibility mode with CUDA 3.2/4.0 and
        # gcc 4.5 compatibility mode with later CUDA versions. This is needed
        # as even with icc used as host compiler, when icc's gcc compatibility
        # mode is higher than the max gcc version officially supported by CUDA,
        # nvcc will freak out.
        if (UNIX AND CMAKE_C_COMPILER_ID MATCHES "Intel" AND
            CUDA_HOST_COMPILER_AUTOSET)
            if (CUDA_VERSION VERSION_LESS "4.1")
                message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.4 for nvcc host compilation")
                set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS};-Xcompiler;-gcc-version=440;")
            else()
                message(STATUS "Setting Intel Compiler compatibity mode to gcc 4.5 for nvcc host compilation")
                set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS};-Xcompiler;-gcc-version=450;")
            endif()
        endif()
        set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS}"
            CACHE STRING "Options for nvcc host compiler (do not edit!)." FORCE)

        mark_as_advanced(CUDA_HOST_COMPILER CUDA_HOST_COMPILER_OPTIONS)
    endif()

    if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "GNU")
        # Some versions of gcc-4.8 and gcc-4.9 produce errors (in particular on OS X)
        # if we do not use -D__STRICT_ANSI__. It is harmless, so we might as well add it for all versions.
        set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS}-D__STRICT_ANSI__;")
    endif()

    # on Linux we need to add -fPIC when building shared gmx libs
    # Note: will add -fPIC for any compiler that supports it as it shouldn't hurt
    if(BUILD_SHARED_LIBS)
        GMX_TEST_CXXFLAG(CXXFLAG_FPIC "-fPIC" _FPIC_NVCC_FLAG)
        if(_FPIC_NVCC_FLAG)
            set(CUDA_HOST_COMPILER_OPTIONS "${CUDA_HOST_COMPILER_OPTIONS}-Xcompiler;${_FPIC_NVCC_FLAG}")
        endif()
    endif()

    # the legacy CUDA kernels have been dropped, warn with CUDA 4.0
    if (CUDA_VERSION VERSION_EQUAL "4.0")
        message(WARNING "The legacy GPU kernels optimized for older CUDA compilers, including the detected version 4.0, have been removed. To avoid performance loss, we strongly recommend upgrading to a newer CUDA toolkit.
        ")
    endif()

    # Set the CUDA GPU architectures to compile for:
    # - with CUDA >v4.2 compute capability 2.0, 2.1 is, but 3.0 is not supported:
    #     => compile sm_20, sm_21 cubin, and compute_20 PTX
    # - with CUDA >=4.2 compute capability <=3.0 is supported:
    #     => compile sm_20, sm_21, sm_30 cubin, and compute_30 PTX
    # - with CUDA 5.0 and later compute capability 3.5 is supported
    #     => compile sm_20, sm_21, sm_30, sm_35 cubin, and compute_35 PTX
    if(CUDA_VERSION VERSION_LESS "4.2")
        set(_CUDA_ARCH_STR "-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_20,code=compute_20")
    elseif(CUDA_VERSION VERSION_LESS "5.0")
        set(_CUDA_ARCH_STR "-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_30,code=sm_30;-gencode;arch=compute_30,code=compute_30")
    else()
        set(_CUDA_ARCH_STR "-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_30,code=sm_30;-gencode;arch=compute_35,code=sm_35;-gencode;arch=compute_35,code=compute_35")
    endif()

    # finally set the damn flags
    set(CUDA_NVCC_FLAGS
        "${_CUDA_ARCH_STR};-use_fast_math;${_HOST_COMPILER_OPTION_STRING}${CUDA_HOST_COMPILER_OPTIONS}"
        CACHE STRING "Compiler flags for nvcc." FORCE)
endif()
