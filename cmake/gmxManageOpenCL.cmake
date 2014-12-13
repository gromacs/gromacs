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

option(GMX_OPENCL_FORCE_CL11_API "Try this if you are having compilations issues with OpenCL enabled" OFF)
option(GMX_OPENCL_HIDE_COMMENT_WARNING "Tell compiler to hide warnings for comments caused by cl_gl_ext.h on Linux" ON)

if(GMX_DOUBLE)
    message(FATAL_ERROR "OpenCL not available in double precision - Yet!")
endif()

#Look for OpenCL
find_package(OpenCL REQUIRED)


#Well the package is REQUIRED so if not found we have stopped already
#But if it was not required consider handling the cases where:
#1) nothing was found (mark something as off and jump out)
#2) OpenCL library was found, but not the headers. Ask the use to install an SDK/Opencl dev package

message(STATUS "OPENCL_INCLUDE_DIRS: " "${OPENCL_INCLUDE_DIRS} ")
message(STATUS "OPENCL_LIBRARIES: " "${OPENCL_LIBRARIES} ")
# detect OpenCL devices in the build host machine
# TO DO: Test the WIN32 branch on Linux
# TO DO: Have just one branch that would work for both Windows and Linux
#if (NOT GMX_OPENCL_DETECTION_DONE)
#	include(gmxDetectGpu)
#	if (WIN32 OR UNIX)
#		gmx_find_OpenCL()
#	else()
#		gmx_detect_OpenCL()
#	endif()
#endif()

#Now configure options
message(STATUS "Setting OpenCL specific options")
    #Where can OpenCL headers be? and with what priority?
    #1: In system
    #2: In paths indicated by environtment variables
    #3: In standard installation paths (e.g. /opt/AMDAPP, /usr/local/cuda etc..
    #4: In Gromacs

if(GMX_OPENCL_FORCE_CL11_API)
    set(OPENCL_DEFINITIONS "-DCL_USE_DEPRECATED_OPENCL_1_1_APIS")
endif(GMX_OPENCL_FORCE_CL11_API)

if(UNIX AND GMX_OPENCL_HIDE_COMMENT_WARNING)
    set(OPENCL_DEFINITIONS ${OPENCL_DEFINITIONS} " -Wno-comment")
endif()

add_definitions(${OPENCL_DEFINITIONS})
include_directories(${OPENCL_INCLUDE_DIRS})

message(STATUS "OpenCL lib: " ${OPENCL_LIBRARIES} ", PATH: " ${OPENCL_INCLUDE_DIRS} ", DEFINITIONS: " ${OPENCL_DEFINITIONS})

macro(gmx_gpu_setup)
    # no OpenMP is no good!
    if(NOT GMX_OPENMP)
        message(WARNING "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading. Without OpenMP a single CPU core can be used with a GPU which is not optimal. Note that with MPI multiple processes can be forced to use a single GPU, but this is typically inefficient. You need to set both C and C++ compilers that support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
    endif()
endmacro()
