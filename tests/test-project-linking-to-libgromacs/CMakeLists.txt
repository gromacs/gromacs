#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2026- The GROMACS Authors
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

cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)

# The GROMACS convention is that these are the version number of the next
# release that is going to be made from this branch.
project(Gromacs-test VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

find_package(Gromacs)

# Set up an executable that links to libgromacs
add_executable(test-gromacs test.cpp)
target_link_libraries(test-gromacs Gromacs::libgromacs)

# Ordinarily the above is enough to test that linking to libgromacs
# suceeds. However, we have had an issue with broken INTERFACE
# properties on imported targets upon which libgromacs depends
# (#5668), but those are usually not a problem. From CMake 4.3, IDEs
# can use CMake's File API to enumerate such targets for their own
# purposes, but we currently only test older versions in GROMACS
# CI. Instead, we can enumerate such targets here and add an explicit
# redundant dependency to expose the problem here.

get_property(_targets DIRECTORY PROPERTY IMPORTED_TARGETS)

foreach(_t IN LISTS _targets)
    get_target_property(_target_type ${_t} TYPE)
    if (${_target_type} STREQUAL "EXECUTABLE")
        # Don't try to link to targets that are executables, it
        # produces a "target not found" error.
    else()
        target_link_libraries(test-gromacs ${_t})
    endif()
endforeach()
