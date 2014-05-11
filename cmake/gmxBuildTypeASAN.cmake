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

# Custom build type "ASAN", to be used to run the
# AddressSanatizer memory checker.

set(_flags "-O1 -g -fsanitize=address -fno-omit-frame-pointer")

foreach(_language C CXX)
    string(REPLACE "X" "+" _human_readable_language ${_language})
    set(CMAKE_${_language}_FLAGS_ASAN ${_flags} CACHE STRING "${_human_readable_language} flags for address sanitizer")
    mark_as_advanced(CMAKE_${_language}_FLAGS_ASAN)
endforeach()

string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
if (APPLE AND _cmake_build_type STREQUAL ASAN) #https://code.google.com/p/address-sanitizer/issues/detail?id=210
   set(BUILD_SHARED_LIBS OFF CACHE BOOL "Disabled for ASAN builds" FORCE)
endif()

