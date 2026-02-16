#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2017- The GROMACS Authors
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

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "Clang is required with GMX_CLANG_CUDA=ON!")
endif()


if(DEFINED GMX_CUDA_CLANG_FLAGS)
    list(APPEND GMX_CUDA_FLAGS ${GMX_CUDA_CLANG_FLAGS})
endif()
# Don't warn about unknown CUDA version; this is developer-facing build, and we hope developers know what they are doing
gmx_add_cuda_flag_if_supported(HAS_NO_UNKNOWN_CUDA_VERSION -Wno-unknown-cuda-version)

# Default flags
gmx_add_cuda_flag_if_supported(HAS_FFAST_MATH -ffast_math)
gmx_add_cuda_flag_if_supported(HAS_FCUDA_FLUSH_DENORMALS_TO_ZERO_FLUSH_DENORMALS_TO_ZERO -fcuda-flush-denormals-to-zero)
