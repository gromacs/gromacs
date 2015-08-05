#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

if(GMX_DOUBLE)
    message(FATAL_ERROR "OpenCL not available in double precision - Yet!")
endif()

# Look for OpenCL
# TODO: FindOpenCL module is available in cmake starting with version 3.1.0.
# A modified version of that module is used here.
# Remove FindOpenCL.cmake file when GROMACS switches to cmake 3.1.0 or higher.
find_package(OpenCL)

if (OPENCL_FOUND)
    if (OPENCL_VERSION_STRING VERSION_LESS REQUIRED_OPENCL_MIN_VERSION)
        message(FATAL_ERROR "OpenCL " "${OPENCL_VERSION_STRING}" " is not supported. OpenCL version " "${REQUIRED_OPENCL_MIN_VERSION}" " or newer is required.")
        return ()
    endif()
else ()
    message(FATAL_ERROR "OpenCL not found.")
    return()
endif()

# Tell compiler to hide warnings for comments caused by cl_gl_ext.h on Linux
if (UNIX)
    set(OPENCL_DEFINITIONS ${OPENCL_DEFINITIONS} " -Wno-comment")
endif()

# Yes Virginia, Darwin kernel version 14.4 corresponds to OS X 10.4.
if(APPLE AND ${CMAKE_SYSTEM_VERSION} VERSION_LESS "14.4")
        message(WARNING "OS X prior to version 10.10.4 produces incorrect AMD OpenCL code at runtime. You will not be able to use AMD GPUs on this host unless you upgrade your operating system.");
endif()

add_definitions(${OPENCL_DEFINITIONS})

include_directories(${OPENCL_INCLUDE_DIRS})

macro(gmx_gpu_setup)
    # no OpenMP is no good!
    if(NOT GMX_OPENMP)
        message(WARNING "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading. Without OpenMP a single CPU core can be used with a GPU which is not optimal. Note that with MPI multiple processes can be forced to use a single GPU, but this is typically inefficient. You need to set both C and C++ compilers that support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
    endif()
endmacro()
