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

# - Try to find TNG
# Once done this will define
#  TNG_FOUND - System has TNG
#  TNG_INCLUDE_DIRS - The TNG include directories
#  TNG_LIBRARIES - The libraries needed to use TNG

find_path(TNG_INCLUDE_DIR tng_io.h PATHS /usr/local/include)

find_library(TNG_LIBRARY NAMES tng_io)

set(TNG_LIBRARIES ${TNG_LIBRARY} )
set(TNG_INCLUDE_DIRS ${TNG_INCLUDE_DIR} )

# if(TNG_INCLUDE_DIR AND TNG_LIBRARY)
#     try_run(TESTTNG TESTTNG_COMPILED ${CMAKE_BINARY_DIR}
#         "${CMAKE_SOURCE_DIR}/cmake/TestTNG.c"
#         CMAKE_FLAGS "-DLINK_LIBRARIES=${TNG_LIBRARY}"
#             "-DINCLUDE_DIRECTORIES=${TNG_INCLUDE_DIR}"
#         RUN_OUTPUT_VARIABLE TNG_VERSION_STRING)
# endif()
#
# handle the QUIETLY and REQUIRED arguments and set TNG_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TNG
                                  REQUIRED_VARS TNG_LIBRARY TNG_INCLUDE_DIR)
#                                   VERSION_VAR TNG_VERSION_STRING)

mark_as_advanced(TNG_INCLUDE_DIR TNG_LIBRARY )

