#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# - Define macro to check if DLOPEN is defined
#
#  GMX_TEST_DLOPEN(VARIABLE)
#
#  VARIABLE will be set if dlopen is present in dlfcn.h
#

macro(GMX_TEST_DLOPEN VARIABLE)
  if(NOT DEFINED ${VARIABLE})
    message(STATUS "Checking for dlopen")

    set(CMAKE_REQUIRED_INCLUDES "dlfcn.h")
    # TODO Make a proper find_package for dlopen to find
    # dlfcn.h. The CMake variable CMAKE_DL_LIBS works magically
    # for the library, however.
    set(CMAKE_REQUIRED_LIBRARIES "dl")
    check_c_source_compiles(
      "#include <dlfcn.h>
int main(void) {
  dlopen(0,0);
  return 0;
}" ${VARIABLE})

    if(${VARIABLE})
      message(STATUS "Checking for dlopen - found")
      set(${VARIABLE} 1 CACHE INTERNAL "Result of test for dlopen" FORCE)
    else()
      message(STATUS "Checking for dlopen - not found")
      set(${VARIABLE} 0 CACHE INTERNAL "Result of test for dlopen" FORCE)
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()
endmacro()
