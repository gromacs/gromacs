#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#
# - Find the Pthreads library
# This module searches for the Pthreads library (including the
# pthreads-win32 port).
#
# This module defines these variables:
#
#  PTHREADS_FOUND
#      True if the Pthreads library was found
#  PTHREADS_LIBRARY
#      The location of the Pthreads library
#  PTHREADS_INCLUDE_DIR
#      The include path of the Pthreads library
#  PTHREADS_DEFINITIONS
#      Preprocessor definitions to define
#
# This module responds to the PTHREADS_EXCEPTION_SCHEME
# variable on Win32 to allow the user to control the
# library linked against.  The Pthreads-win32 port
# provides the ability to link against a version of the
# library with exception handling.  IT IS NOT RECOMMENDED
# THAT YOU USE THIS because most POSIX thread implementations
# do not support stack unwinding.
#
#  PTHREADS_EXCEPTION_SCHEME
#      C  = no exceptions (default)
#         (NOTE: This is the default scheme on most POSIX thread
#          implementations and what you should probably be using)
#      CE = C++ Exception Handling
#      SE = Structure Exception Handling (MSVC only)
#

#
# Define a default exception scheme to link against
# and validate user choice.
#
IF(PTHREADS_INCLUDE_DIR)
  # Already in cache, be silent
  SET(PTHREADS_FIND_QUIETLY TRUE)
ENDIF(PTHREADS_INCLUDE_DIR)


IF(NOT DEFINED PTHREADS_EXCEPTION_SCHEME)
    # Assign default if needed
    SET(PTHREADS_EXCEPTION_SCHEME "C")
ELSE(NOT DEFINED PTHREADS_EXCEPTION_SCHEME)
    # Validate
    IF(NOT PTHREADS_EXCEPTION_SCHEME STREQUAL "C" AND
       NOT PTHREADS_EXCEPTION_SCHEME STREQUAL "CE" AND
       NOT PTHREADS_EXCEPTION_SCHEME STREQUAL "SE")

    MESSAGE(FATAL_ERROR "See documentation for FindPthreads.cmake, only C, CE, and SE modes are allowed")

    ENDIF(NOT PTHREADS_EXCEPTION_SCHEME STREQUAL "C" AND
          NOT PTHREADS_EXCEPTION_SCHEME STREQUAL "CE" AND
          NOT PTHREADS_EXCEPTION_SCHEME STREQUAL "SE")

     IF(NOT MSVC AND PTHREADS_EXCEPTION_SCHEME STREQUAL "SE")
         MESSAGE(FATAL_ERROR "Structured Exception Handling is only allowed for MSVC")
     ENDIF(NOT MSVC AND PTHREADS_EXCEPTION_SCHEME STREQUAL "SE")

ENDIF(NOT DEFINED PTHREADS_EXCEPTION_SCHEME)

#
# Find the header file
#
FIND_PATH(PTHREADS_INCLUDE_DIR pthread.h)

#
# Find the library
#
SET(names)
IF(MSVC)
    SET(names
            pthreadV${PTHREADS_EXCEPTION_SCHEME}2
            pthread
    )
ELSEIF(MINGW)
    SET(names
            pthreadG${PTHREADS_EXCEPTION_SCHEME}2
            pthread
    )
ELSE(MSVC) # Unix / Cygwin / Apple
    SET(names pthread)
ENDIF(MSVC)
    
FIND_LIBRARY(PTHREADS_LIBRARY ${names}
    DOC "The Portable Threads Library")

IF(PTHREADS_INCLUDE_DIR AND PTHREADS_LIBRARY)
    SET(PTHREADS_FOUND true)
    SET(PTHREADS_DEFINITIONS -DHAVE_PTHREAD_H)
    SET(PTHREADS_INCLUDE_DIRS ${PTHREADS_INCLUDE_DIR})
    SET(PTHREADS_LIBRARIES    ${PTHREADS_LIBRARY})
ENDIF(PTHREADS_INCLUDE_DIR AND PTHREADS_LIBRARY)

IF(PTHREADS_FOUND)
    IF(NOT PTHREADS_FIND_QUIETLY)
        MESSAGE(STATUS "Found Pthreads: ${PTHREADS_LIBRARY}")
    ENDIF(NOT PTHREADS_FIND_QUIETLY)
ELSE(PTHREADS_FOUND) 
    IF(PTHREADS_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find the Pthreads Library")
    ENDIF(PTHREADS_FIND_REQUIRED)
ENDIF(PTHREADS_FOUND)
