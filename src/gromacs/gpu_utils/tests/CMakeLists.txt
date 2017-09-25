#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

# Arrange to compile files with test infrastructure as normal C++, and
# other files that sometimes have platform-specific symbols using the
# required infrastructure to compile and link, without needing three
# different files and strategies for the different GPU implementation
# flavours.

# Always compiled as plain C++
file(GLOB SOURCES_FROM_CXX
    hostallocator.cpp
    )

# Compiled with nvcc if necessary
file(GLOB GPU_SOURCES_FROM_CXX
    hostallocator-impl.cpp
    )

if(GMX_USE_CUDA)
    # Compile source files with .cpp extension with nvcc using -x cu
    foreach(_cuda_file ${GPU_SOURCES_FROM_CXX})
        set_source_files_properties(${_cuda_file} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    endforeach()
    # TODO Making a separate library is heavy handed, but nothing else
    # seems to work. Also don't use a hyphen in its name, because nvcc
    # can't cope with that.
    cuda_add_library(libgpu_utilstest ${GPU_SOURCES_FROM_CXX})
else()
    # Do normal compilation of files that might have been necessary to
    # compile with nvcc.
    list(APPEND SOURCES_FROM_CXX ${GPU_SOURCES_FROM_CXX})
endif()

gmx_add_unit_test(GpuUtilsUnitTests gpu_utils-test
    # Infrastructure
    gputest.cpp
    # Tests of code
    ${SOURCES_FROM_CXX}
    )

if(GMX_USE_CUDA)
    target_link_libraries(gpu_utils-test libgpu_utilstest)
endif()
if(GMX_USE_OPENCL)
    target_link_libraries(gpu_utils-test ${OPENCL_LIBRARIES})
endif()
