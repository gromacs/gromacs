#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
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


gmx_check_if_changed(CUDA_HOST_COMPILER_CHANGED CMAKE_CUDA_HOST_COMPILER)
# glibc 2.23 changed string.h in a way that breaks CUDA compilation in
# many projects, but which has a trivial workaround. It would be nicer
# to compile with nvcc and see that the workaround is necessary and
# effective, but it is unclear how to do that. Also, grepping in the
# glibc source shows that _FORCE_INLINES is only used in this string.h
# feature and performance of memcpy variants is unimportant for CUDA
# code in GROMACS. So this workaround is good enough to keep problems
# away from users installing GROMACS. See Issue #1982.
if(CUDA_HOST_COMPILER_CHANGED)
    try_compile(IS_GLIBC_2_23_OR_HIGHER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/TestGlibcVersion.cpp)
    if(IS_GLIBC_2_23_OR_HIGHER)
        message(STATUS "Adding work-around for issue compiling CUDA code with glibc 2.23 string.h")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D_FORCE_INLINES")
    endif()
endif()

# pass the GMXC_CXXFLAGS to nvcc as host compiler options.
# This is used to suppress warnings for Clang for now.
# Ideally should be enabled for all host compilers once
# all the warnings are resolved in the test cases which are
# compiled with nvcc.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(GMXC_CXXFLAGS)
        # convert the string to list
        string(REPLACE " " ";" GMXC_CXXFLAGS_LIST_ "${GMXC_CXXFLAGS}")
        # Loop over each item in the list
        foreach(item IN LISTS GMXC_CXXFLAGS_LIST_)
            list(APPEND CUDA_HOST_COMPILER_OPTIONS -Xcompiler="${item}")
        endforeach()
    endif()
endif()

# assemble the CUDA flags
gmx_add_cuda_flag_if_supported(NVCC_HAS_USE_FAST_MATH -use_fast_math)
gmx_add_cuda_flag_if_supported(NVCC_HAS_STATIC_GLOBAL_TEMPLATE_STUB_FALSE -static-global-template-stub=false)
# Add warnings
gmx_add_cuda_flag_if_supported(NVCC_HAS_PTXAS_WARN_DOUBLE_USAGE -Xptxas=-warn-double-usage)
gmx_add_cuda_flag_if_supported(NVCC_HAS_PTXAS_WERROR -Xptxas=-Werror)

# Disable cudafe warnings with nvc++ as a host compiler - warning #177-D
gmx_add_cuda_flag_if_supported(NVCC_HAS_DIAG_SUPPRESS_177 -diag-suppress=177)

gmx_check_if_changed(_cuda_nvcc_executable_or_flags_changed CMAKE_CUDA_COMPILER GMX_CUDA_FLAGS)
