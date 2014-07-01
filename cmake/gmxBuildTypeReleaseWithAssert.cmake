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

# Custom build type "ReleaseWithAssert".
string(REGEX REPLACE "[/-][dD][^/-]*NDEBUG" "" CMAKE_C_FLAGS_RELWITHASSERT_INIT "${CMAKE_C_FLAGS_RELEASE_INIT}")
string(REGEX REPLACE "[/-][dD][^/-]*NDEBUG" "" CMAKE_CXX_FLAGS_RELWITHASSERT_INIT "${CMAKE_CXX_FLAGS_RELEASE_INIT}")
set( CMAKE_C_FLAGS_RELWITHASSERT "${CMAKE_C_FLAGS_RELWITHASSERT_INIT}" CACHE STRING "C flags for release with assert builds.")
set( CMAKE_CXX_FLAGS_RELWITHASSERT "${CMAKE_CXX_FLAGS_RELWITHASSERT_INIT}" CACHE STRING "C++ flags for release with assert builds.")
set( CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Linker flags for release with assert builds.")
mark_as_advanced( CMAKE_CXX_FLAGS_RELWITHASSERT CMAKE_C_FLAGS_RELWITHASSERT)

