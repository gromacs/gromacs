#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016,2017,2020, by the GROMACS development team, led by
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

gmx_add_libgromacs_sources(
    cpuinfo.cpp
    detecthardware.cpp
    device_management_common.cpp
    hardwaretopology.cpp
    printhardware.cpp
    identifyavx512fmaunits.cpp
    )

if(GMX_GPU_OPENCL)
    gmx_add_libgromacs_sources(
        device_management_ocl.cpp
        )
elseif(GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
        device_management.cu
        )
    _gmx_add_files_to_property(CUDA_SOURCES
        # Must add these files so they can include device_information.h
        device_management_common.cpp
        detecthardware.cpp
	)
elseif(GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        device_management_sycl.cpp
    )
    _gmx_add_files_to_property(SYCL_SOURCES
        device_management_sycl.cpp
        # Must add these files so they can include device_information.h
        device_management.cpp
        device_management_common.cpp
        detecthardware.cpp
	)
else()
    gmx_add_libgromacs_sources(
        device_management.cpp
    )
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
