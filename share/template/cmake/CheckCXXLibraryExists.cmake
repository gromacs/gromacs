#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

# - Check if the function exists in a C++ library.
# CHECK_CXX_LIBRARY_EXISTS (LIBRARY FUNCTION LOCATION VARIABLE)
#
#  LIBRARY  - the name of the library you are looking for
#  FUNCTION - the name of the function
#  LOCATION - location where the library should be found
#  VARIABLE - variable to store the result
#
# The following variables may be set before calling this macro to
# modify the way the check is run:
#
#  CMAKE_REQUIRED_FLAGS = string of compile command line flags
#  CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#  CMAKE_REQUIRED_LIBRARIES = list of libraries to link

#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
#
# CMake - Cross Platform Makefile Generator
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.

# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# ------------------------------------------------------------------------------

# The above copyright and license notice applies to distributions of
# CMake in source and binary form.  Some source files contain additional
# notices of original copyright by their contributors; see each source
# for details.  Third-party software packages supplied with CMake under
# compatible licenses provide their own copyright notices documented in
# corresponding subdirectories.

# ------------------------------------------------------------------------------

# CMake was initially developed by Kitware with the following sponsorship:

#  * National Library of Medicine at the National Institutes of Health
#    as part of the Insight Segmentation and Registration Toolkit (ITK).

#  * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#    Visualization Initiative.

#  * National Alliance for Medical Image Computing (NAMIC) is funded by the
#    National Institutes of Health through the NIH Roadmap for Medical Research,
#    Grant U54 EB005149.

#  * Kitware, Inc.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# Imported by Mark Abraham from CheckLibraryExists.cmake in cmake
# 2.8.8, and modified for use with FindGROMACS, to detect a function
# in a C++ library via try_compile of a C++ program rather than a C
# program. In particular, the use of LINK_LIBRARIES below changed in
# cmake 2.8.11, so this could be updated, when the GROMACS (template)
# minimum required version is increased.

INCLUDE(CMakeExpandImportedTargets)

MACRO(CHECK_CXX_LIBRARY_EXISTS LIBRARY FUNCTION_DECLARATION FUNCTION_CALL LOCATION VARIABLE)
  IF("${VARIABLE}" MATCHES "^${VARIABLE}$")
    SET(MACRO_CHECK_LIBRARY_EXISTS_DEFINITION
      "-DCHECK_CXX_FUNCTION_EXISTS_DECLARATION='${FUNCTION_DECLARATION}' -DCHECK_CXX_FUNCTION_EXISTS_CALL='${FUNCTION_CALL}' ${CMAKE_REQUIRED_FLAGS}")
    MESSAGE(STATUS "Looking for ${FUNCTION_DECLARATION} in ${LIBRARY}")
    SET(CHECK_LIBRARY_EXISTS_LIBRARIES ${LIBRARY})
    IF(CMAKE_REQUIRED_LIBRARIES)
      # this one translates potentially used imported library targets to their files on disk
      CMAKE_EXPAND_IMPORTED_TARGETS(_ADJUSTED_CMAKE_REQUIRED_LIBRARIES  LIBRARIES  ${CMAKE_REQUIRED_LIBRARIES} CONFIGURATION "${CMAKE_TRY_COMPILE_CONFIGURATION}")
      SET(CHECK_LIBRARY_EXISTS_LIBRARIES
        ${CHECK_LIBRARY_EXISTS_LIBRARIES} ${_ADJUSTED_CMAKE_REQUIRED_LIBRARIES})
    ENDIF(CMAKE_REQUIRED_LIBRARIES)
    TRY_COMPILE(${VARIABLE}
      ${CMAKE_BINARY_DIR}
      ${CMAKE_CURRENT_LIST_DIR}/CheckCXXFunctionExists.cpp
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      CMAKE_FLAGS
      -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_LIBRARY_EXISTS_DEFINITION}
      -DLINK_DIRECTORIES:STRING=${LOCATION}
      "-DLINK_LIBRARIES:STRING=${CHECK_LIBRARY_EXISTS_LIBRARIES}"
      OUTPUT_VARIABLE OUTPUT)

    IF(${VARIABLE})
      MESSAGE(STATUS "Looking for ${FUNCTION_DECLARATION} in ${LIBRARY} - found")
      SET(${VARIABLE} 1 CACHE INTERNAL "Have library ${LIBRARY}")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Determining if the function ${FUNCTION_DECLARATION} exists in the ${LIBRARY} "
        "passed with the following output:\n"
        "${OUTPUT}\n\n")
    ELSE(${VARIABLE})
      MESSAGE(STATUS "Looking for ${FUNCTION_DECLARATION} in ${LIBRARY} - not found")
      SET(${VARIABLE} "" CACHE INTERNAL "Have library ${LIBRARY}")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Determining if the function ${FUNCTION_DECLARATION} exists in the ${LIBRARY} "
        "failed with the following output:\n"
        "${OUTPUT}\n\n")
    ENDIF(${VARIABLE})
  ENDIF("${VARIABLE}" MATCHES "^${VARIABLE}$")
ENDMACRO(CHECK_CXX_LIBRARY_EXISTS)
