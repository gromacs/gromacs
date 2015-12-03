#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2015, by the GROMACS development team, led by
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

# - Define macro to check if DLOPEN is defined
#
#  GMX_TEST_DLOPEN(VARIABLE)
#
#  VARIABLE will be set if dlopen is present in dlfcn.h
#

MACRO(GMX_TEST_DLOPEN VARIABLE)
  IF(NOT DEFINED ${VARIABLE})
    MESSAGE(STATUS "Checking for dlopen")

    set(CMAKE_REQUIRED_INCLUDES "dlfcn.h")
    set(CMAKE_REQUIRED_LIBRARIES "dl")
    check_c_source_compiles(
      "#include <dlfcn.h>
int main(void) {
  dlopen(0,0);
  return 0;
}" ${VARIABLE})

    IF(${VARIABLE})
      MESSAGE(STATUS "Checking for dlopen - found")
      set(${VARIABLE} 1 CACHE INTERNAL "Result of test for dlopen" FORCE)
    ELSE()
      MESSAGE(STATUS "Checking for dlopen - not found")
      set(${VARIABLE} 0 CACHE INTERNAL "Result of test for dlopen" FORCE)
    ENDIF()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  ENDIF()
ENDMACRO()
