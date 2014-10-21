#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2014, by the GROMACS development team, led by
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

# Defines macros to detect various functions related to dynamic linker
# interaction
#
# The following macros are provided:
#   gmx_test_dlopen(VARIABLE)
#     Sets VARIABLE to 1 or 0 based on whether dlopen is available.
#   gmx_test_dladdr(VARIABLE)
#     Sets VARIABLE to 1 or 0 based on whether dladdr is available.

include(CheckCSourceCompiles)

function (gmx_test_dlopen VARIABLE)
    # The variable is created in the cache by check_c_source_compiles()
    if (NOT DEFINED ${VARIABLE})
        set(CMAKE_REQUIRED_DEFINITIONS "")
        set(CMAKE_REQUIRED_FLAGS "")
        set(CMAKE_REQUIRED_INCLUDES "")
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})
        check_c_source_compiles("
#include <dlfcn.h>
int main(void) { dlopen(0,0); }" ${VARIABLE})
    endif()
endfunction()

function (gmx_test_dladdr VARIABLE)
    # The variable is created in the cache by check_c_source_compiles()
    if (NOT DEFINED ${VARIABLE})
        set(CMAKE_REQUIRED_DEFINITIONS "")
        set(CMAKE_REQUIRED_FLAGS "")
        set(CMAKE_REQUIRED_INCLUDES "")
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})
        check_c_source_compiles("
#include <dlfcn.h>
int main(void) {
  Dl_info info;
  dladdr(0,&info);
  return (info.dli_fname != 0) ? 0 : 1;
}" ${VARIABLE})
    endif()
endfunction()
