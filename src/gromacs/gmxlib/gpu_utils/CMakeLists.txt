#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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

# (slightly sloppy) OS definitions required by memtestG80
set(_os_def)
if(UNIX)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(_os_def "-DOSX")
    else() # everything that's UNIX & UNIX-like except OS X
        set(_os_def "-DLINUX")
    endif()
else()
    if(WIN32)
        set(_os_def "-DWINDOWS")
    else()
        message(FATAL_ERROR " Could not detect OS required for memtestG80.")
    endif()
endif()

CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)        
file(GLOB GPU_UTILS_SOURCES *.cu)
CUDA_ADD_LIBRARY(gpu_utils STATIC ${GPU_UTILS_SOURCES}
                 OPTIONS ${_os_def}
                 RELWITHDEBINFO -g
                 DEBUG -g -D_DEBUG_=1 )

CUDA_BUILD_CLEAN_TARGET()
