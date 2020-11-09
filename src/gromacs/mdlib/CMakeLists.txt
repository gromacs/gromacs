#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2012,2013,2014,2015 by the GROMACS development team.
# Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

file(GLOB MDLIB_SOURCES *.cpp)
# To avoid listing all the necessary files manually, we will remove SYCL-specfific one here:
list(REMOVE_ITEM MDLIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/leapfrog_gpu_sycl.cpp)

set(MDLIB_SOURCES ${MDLIB_SOURCES} PARENT_SCOPE)
if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
if(GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
       leapfrog_gpu.cu
       lincs_gpu.cu
       settle_gpu.cu
       update_constrain_gpu_impl.cu
       gpuforcereduction_impl.cu
       )
endif()
if(GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        leapfrog_gpu_sycl.cpp
    )
    _gmx_add_files_to_property(SYCL_SOURCES
        leapfrog_gpu_sycl.cpp
    )
endif()
