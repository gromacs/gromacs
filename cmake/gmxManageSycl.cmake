#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
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

set(GMX_GPU_SYCL ON)

if(GMX_DOUBLE)
    message(FATAL_ERROR "SYCL acceleration is not available in double precision")
endif()

set(_sycl_has_valid_fft FALSE)

if(GMX_SYCL STREQUAL "ACPP")
    set(GMX_SYCL_DPCPP OFF)
    set(GMX_SYCL_ACPP ON)
    include(gmxManageSyclAdaptiveCpp)
elseif(GMX_SYCL STREQUAL "DPCPP")
    set(GMX_SYCL_DPCPP ON)
    set(GMX_SYCL_ACPP OFF)
    include(gmxManageSyclOneApi)
else()
    message(FATAL_ERROR "Unsupported value for GMX_SYCL: ${GMX_SYCL}. Please set either \"ACPP\" or \"DPCPP\"")
endif()

if (GMX_GPU_FFT_CUFFT AND GMX_USE_HEFFTE)
    set(_sycl_has_valid_fft TRUE)
    if (NOT DEFINED ENV{GITLAB_CI}) # Don't warn in CI builds
        message(WARNING "SYCL build with HeFFTe and cuFFT should only ever be used for testing")
    endif()
endif()

if(NOT ${_sycl_has_valid_fft} AND NOT GMX_GPU_FFT_LIBRARY STREQUAL "NONE")
    set(_hint "")
    if (GMX_GPU_FFT_CUFFT OR GMX_GPU_FFT_CLFFT)
        set(_hint " It is not supported with SYCL.")
    elseif (GMX_SYCL_ACPP AND GMX_GPU_FFT_MKL)
        set(_hint " MKL is only supported with Intel DPC++ compiler, not with AdaptiveCpp/hipSYCL")
    endif()
    message(FATAL_ERROR "The selected GPU FFT library ${GMX_GPU_FFT_LIBRARY} is not compatible.${_hint}")
endif()

if(NOT ${_sycl_has_valid_fft} AND NOT DEFINED ENV{GITLAB_CI}) # Don't warn in CI builds
    message(WARNING "Building SYCL version without GPU FFT library.  Will not be able to perform FFTs on a GPU, which is not good for performance.")
endif()
