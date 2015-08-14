#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#.rst:
# FindNVML
# --------
#
# Find the NVIDIA Management Library (NVML) includes and library. NVML documentation
# is available at: http://docs.nvidia.com/deploy/nvml-api/index.html 
#
# NVML is part of the GPU Deployment Kit (GDK) and GPU_DEPLOYMENT_KIT_ROOT_DIR can
# be specified if the GPU Deployment Kit is not installed in a default location.
#
# FindNVML defines the following variables: 
#
#   NVML_INCLUDE_DIR, where to find nvml.h, etc.
#   NVML_LIBRARY, the libraries needed to use NVML.
#   NVML_FOUND, If false, do not try to use NVML.
#

#   Jiri Kraus, NVIDIA Corp (nvidia.com - jkraus)
#
#   Copyright (c) 2008 - 2014 NVIDIA Corporation.  All rights reserved.
#
#   This code is licensed under the MIT License.  See the FindNVML.cmake script
#   for the text of the license.

# The MIT License
#
# License for the specific language governing rights and limitations under
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
###############################################################################

if( CMAKE_SYSTEM_NAME STREQUAL "Windows"  )
  set( NVML_LIB_PATHS "C:/Program Files/NVIDIA Corporation/GDK/nvml/lib" )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
    list(APPEND NVML_LIB_PATHS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/nvml/lib")
  endif()
  set(NVML_NAMES nvml)
  
  set( NVML_INC_PATHS "C:/Program Files/NVIDIA Corporation/GDK/nvml/include" )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
    list(APPEND NVML_INC_PATHS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/nvml/include")
  endif()
else()
  # The Linux installer for the GPU Deployment Kit adds a "usr"
  # suffix to a custom path if one is used, so a user could
  # reasonably set GPU_DEPLOYMENT_KIT_ROOT_DIR to the value they
  # passed to the installer, or the root where they later found the
  # kit to be installed. Below, we cater for both possibilities.
  set( NVML_LIB_PATHS /usr/lib64 )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
      list(APPEND NVML_LIB_PATHS
          "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/src/gdk/nvml/lib"
          "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/usr/src/gdk/nvml/lib"
          )
  endif()
  set(NVML_NAMES nvidia-ml)
  
  set( NVML_INC_PATHS /usr/include/nvidia/gdk/ /usr/include )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
      list(APPEND NVML_INC_PATHS
          "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/include/nvidia/gdk"
          "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/usr/include/nvidia/gdk"
          )
  endif()
endif()

find_library(NVML_LIBRARY NAMES ${NVML_NAMES} PATHS ${NVML_LIB_PATHS} )

find_path(NVML_INCLUDE_DIR nvml.h PATHS ${NVML_INC_PATHS})

# handle the QUIETLY and REQUIRED arguments and set NVML_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NVML DEFAULT_MSG NVML_LIBRARY NVML_INCLUDE_DIR)

mark_as_advanced(NVML_LIBRARY NVML_INCLUDE_DIR)
