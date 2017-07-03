#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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

# Set up files that are portably compiled from a single cpp file,
# regardless of what flavour of GPU build this is.
file(GLOB GPU_SOURCES_FROM_CXX
    hostallocator.cpp
    pagelockedmemory.cpp
    )

# Compile files that are non-portable

if(GMX_USE_OPENCL)
    gmx_add_libgromacs_sources(gpu_utils_ocl.cpp ocl_compiler.cpp ocl_caching.cpp oclutils.cpp ${GPU_SOURCES_FROM_CXX})
elseif(GMX_USE_CUDA)
    gmx_add_libgromacs_sources(cudautils.cu gpu_utils.cu pmalloc_cuda.cu ${GPU_SOURCES_FROM_CXX})
    # For CUDA, we later need to set a cmake property on each such
    # file, so nvcc knows that it should compile it with nvcc -x cu,
    # so we pass the full list back to the parent CMakeLists.txt.
    list(APPEND CUDA_SOURCES_FROM_CXX ${GPU_SOURCES_FROM_CXX})
    set(CUDA_SOURCES_FROM_CXX ${CUDA_SOURCES_FROM_CXX} PARENT_SCOPE)
else()
    gmx_add_libgromacs_sources(gpu_utils.cpp ${GPU_SOURCES_FROM_CXX})
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
