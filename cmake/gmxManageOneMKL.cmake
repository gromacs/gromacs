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
# 
# Copyright 2024- Codeplay Software Ltd.

# Manage OneMKL, GPU FFT library used with SYCL.

function(gmx_manage_onemkl)
    # Find quietly the second time.
    if (oneMKL_FIND_QUIETLY_AFTER_FIRST_RUN)
        set(oneMKL_FIND_QUIETLY TRUE)
    endif()
    find_package(oneMKL REQUIRED)
    set(oneMKL_FIND_QUIETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to find oneMKL")

    set(BACKEND_COUNT 0)
    if("${SYCL_CXX_FLAGS_EXTRA}" MATCHES "fsycl-targets=.*(nvptx64|nvidia_gpu)")
        if(NOT TARGET MKL::onemkl_dft_cufft)
            message(WARNING "GROMACS SYCL is targetting NVIDIA GPU, but oneMKL interface library was not built with the cuFFT backend.")
        else()
            MATH(EXPR BACKEND_COUNT "${BACKEND_COUNT}+1")
        endif()
    elseif("${SYCL_CXX_FLAGS_EXTRA}" MATCHES "fsycl-targets=.*(amdgcn|amd_gpu)")
        if(NOT TARGET MKL::onemkl_dft_rocfft)
            message(WARNING "GROMACS SYCL is targetting AMD GPU, but oneMKL interface library was not built with the rocFFT backend.")
        else()
            MATH(EXPR BACKEND_COUNT "${BACKEND_COUNT}+1")
        endif()
    else()
        if(NOT TARGET MKL::onemkl_dft_mklgpu)
            message(WARNING "No oneMKL interface library DFT backend is compatible with the chosen fsycl-targets.")
        else()
            MATH(EXPR BACKEND_COUNT "${BACKEND_COUNT}+1")
        endif()
    endif()

    if(${BACKEND_COUNT} EQUAL 0)
        message(FATAL_ERROR "oneMKL interface library was not compatible with -fsycl-targets")
    endif()
    if(${BACKEND_COUNT} GREATER 1)
        message(FATAL_ERROR "GROMACS only supports linking against a single oneMKL backend. Too many compatible fsycl-targets.")
    endif()
endfunction()
