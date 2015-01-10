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

#.rst:
# FindOPENCL
# ----------
#
# Try to find OPENCL
#
# Once done this will define::
#
#   OPENCL_FOUND          - True if OPENCL was found
#   OPENCL_INCLUDE_DIRS   - include directories for OPENCL
#   OPENCL_LIBRARIES      - link against this library to use OPENCL
#   OPENCL_VERSION_STRING - Highest supported OPENCL version (eg. 1.2)
#   OPENCL_VERSION_MAJOR  - The major version of the OPENCL implementation
#   OPENCL_VERSION_MINOR  - The minor version of the OPENCL implementation
#
# The module will also define two cache variables::
#
#   OPENCL_INCLUDE_DIR    - the OPENCL include directory
#   OPENCL_LIBRARY        - the path to the OPENCL library
#
# This is a modified version of FindOpenCL.cmake from cmake v3.1.0
# (see comments at the end of the file).
# The following changes have been made:
#     1. OpenCL is written in all caps (OPENCL)
#     2. The following block has been modified:
#include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
#find_package_handle_standard_args(
#  OpenCL
#  FOUND_VAR OpenCL_FOUND
#  REQUIRED_VARS OpenCL_LIBRARY OpenCL_INCLUDE_DIR
#  VERSION_VAR OpenCL_VERSION_STRING)
#     has been replaced by:
#include(FindPackageHandleStandardArgs)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENCL
#  REQUIRED_VARS OPENCL_LIBRARY OPENCL_INCLUDE_DIR
#  VERSION_VAR OPENCL_VERSION_STRING)
#     3. The following block has been modified:
#  find_library(OPENCL_LIBRARY
#   NAMES OPENCL)
#     has been replaced by:
#  find_library(OPENCL_LIBRARY
#   NAMES OpenCL)

function(_FIND_OPENCL_VERSION)
  include(CheckSymbolExists)
  include(CMakePushCheckState)
  set(CMAKE_REQUIRED_QUIET ${OPENCL_FIND_QUIETLY})

  CMAKE_PUSH_CHECK_STATE()
  foreach(VERSION "2_0" "1_2" "1_1" "1_0")
    set(CMAKE_REQUIRED_INCLUDES "${OPENCL_INCLUDE_DIR}")

    if(APPLE)
      CHECK_SYMBOL_EXISTS(
        CL_VERSION_${VERSION}
        "OpenCL/cl.h"
        OPENCL_VERSION_${VERSION})
    else()
      CHECK_SYMBOL_EXISTS(
        CL_VERSION_${VERSION}
        "CL/cl.h"
        OPENCL_VERSION_${VERSION})
    endif()

    if(OPENCL_VERSION_${VERSION})
      string(REPLACE "_" "." VERSION "${VERSION}")
      set(OPENCL_VERSION_STRING ${VERSION} PARENT_SCOPE)
      string(REGEX MATCHALL "[0-9]+" version_components "${VERSION}")
      list(GET version_components 0 major_version)
      list(GET version_components 1 minor_version)
      set(OPENCL_VERSION_MAJOR ${major_version} PARENT_SCOPE)
      set(OPENCL_VERSION_MINOR ${minor_version} PARENT_SCOPE)
      break()
    endif()
  endforeach()
  CMAKE_POP_CHECK_STATE()
endfunction()

find_path(OPENCL_INCLUDE_DIR
  NAMES
    CL/cl.h OpenCL/cl.h
  PATHS
    ENV "PROGRAMFILES(X86)"
    ENV AMDAPPSDKROOT
    ENV INTELOCLSDKROOT
    ENV NVSDKCOMPUTE_ROOT
    ENV CUDA_PATH
    ENV CUDA_HOME
    ENV ATISTREAMSDKROOT
  PATH_SUFFIXES
    include
    OPENCL/common/inc
    "AMD APP/include")

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    find_library(OPENCL_LIBRARY
        NAMES OPENCL OpenCL
        PATHS
        ENV "PROGRAMFILES(X86)"
        ENV AMDAPPSDKROOT
        ENV INTELOCLSDKROOT
        ENV CUDA_PATH
        ENV CUDA_HOME
        ENV NVSDKCOMPUTE_ROOT
        ENV ATISTREAMSDKROOT
        PATH_SUFFIXES
        "AMD APP/lib/x86"
        lib/x86
        lib/Win32
        lib
        OPENCL/common/lib/Win32)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
    find_library(OPENCL_LIBRARY
        NAMES OPENCL OpenCL
        PATHS
        ENV "PROGRAMFILES(X86)"
        ENV AMDAPPSDKROOT
        ENV INTELOCLSDKROOT
        ENV CUDA_PATH
        ENV CUDA_HOME
        ENV NVSDKCOMPUTE_ROOT
        ENV ATISTREAMSDKROOT
        PATH_SUFFIXES
        "AMD APP/lib/x86_64"
        lib/x86_64
        lib/x64
        lib64
        OPENCL/common/lib/x64)
endif()

_FIND_OPENCL_VERSION()

set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
set(OPENCL_INCLUDE_DIRS ${OPENCL_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENCL
  REQUIRED_VARS OPENCL_LIBRARY OPENCL_INCLUDE_DIR
  VERSION_VAR OPENCL_VERSION_STRING)

mark_as_advanced(
  OPENCL_INCLUDE_DIR
  OPENCL_LIBRARY)
 
#=============================================================================
# Copyright 2014 Matthaeus G. Chajdas
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# cmake 3.1.0 Copyright.txt file content is attached below:
#
#CMake - Cross Platform Makefile Generator
#Copyright 2000-2014 Kitware, Inc.
#Copyright 2000-2011 Insight Software Consortium
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#
#* Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
#* Neither the names of Kitware, Inc., the Insight Software Consortium,
#  nor the names of their contributors may be used to endorse or promote
#  products derived from this software without specific prior written
#  permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#------------------------------------------------------------------------------
#
#The above copyright and license notice applies to distributions of
#CMake in source and binary form.  Some source files contain additional
#notices of original copyright by their contributors; see each source
#for details.  Third-party software packages supplied with CMake under
#compatible licenses provide their own copyright notices documented in
#corresponding subdirectories.
#
#------------------------------------------------------------------------------
#
#CMake was initially developed by Kitware with the following sponsorship:
#
# * National Library of Medicine at the National Institutes of Health
#   as part of the Insight Segmentation and Registration Toolkit (ITK).
#
# * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#   Visualization Initiative.
#
# * National Alliance for Medical Image Computing (NAMIC) is funded by the
#   National Institutes of Health through the NIH Roadmap for Medical Research,
#   Grant U54 EB005149.
#
# * Kitware, Inc.
