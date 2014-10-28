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

include(Platform/Linux-Intel-Cray)
__linux_compiler_intel(C)
__linux_compiler_intel(CXX)
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")

# Set up cache variables that should always be this way (no need to
# use FORCE, the user might know what they are doing).
set(CMAKE_SKIP_RPATH ON CACHE BOOL "Can't use RPATH with Cray")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Building statically on Cray")
set(GMX_PREFER_STATIC_LIBS ON CACHE BOOL "Building statically on Cray")
set(GMX_LINK_STATIC_BINARIES ON CACHE BOOL "Avoid linking shared versions of e.g. system libraries")
