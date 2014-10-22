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

#TODO: Should GPU_DEPLOYMENT_KIT_ROOT_DIR override the default locations?
if( CMAKE_SYSTEM_NAME STREQUAL "Windows"  )
  set( NVML_LIB_PATHS "C:\Program Files\NVIDIA Corporation\GDK" )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
    list(APPEND NVML_LIB_PATHS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}\nvml\lib")
  endif()
  set(NVML_NAMES nvml)
  
  set( NVML_INC_PATHS "C:\Program Files\NVIDIA Corporation\GDK\nvml\include" )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
    list(APPEND NVML_INC_PATHS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}\nvml\include")
  endif()
else()
  set( NVML_LIB_PATHS /usr/lib64 )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
    list(APPEND NVML_LIB_PATHS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/src/gdk/nvml/lib")
  endif()
  set(NVML_NAMES nvidia-ml)
  
  set( NVML_INC_PATHS /usr/include/nvidia/gdk/ /usr/include )
  if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
    list(APPEND NVML_INC_PATHS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/include/nvidia/gdk")
  endif()
endif()

find_library(NVML_LIBRARY NAMES ${NVML_NAMES} PATHS ${NVML_LIB_PATHS} )

find_path(NVML_INCLUDE_DIR nvml.h PATHS ${NVML_INC_PATHS})

# handle the QUIETLY and REQUIRED arguments and set NVML_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NVML DEFAULT_MSG NVML_LIBRARY NVML_INCLUDE_DIR)

mark_as_advanced(NVML_LIBRARY NVML_INCLUDE_DIR)
