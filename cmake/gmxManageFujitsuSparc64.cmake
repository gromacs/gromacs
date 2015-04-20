#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# Managing configuration for Fujitsu PrimeHPC Sparc64
# For now this is mainly used for K computer.
message(STATUS "Configuring for Fujitsu Sparc64")

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Use static linking by default on Fujitsu Sparc64" FORCE)
set(GMX_GPU OFF CACHE BOOL "Cannot do GPU acceleration on Fujitsu Sparc64" FORCE)

set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Use native 1.0/sqrt(x) on Fujitsu Sparc64" FORCE)
set(GMX_X11 OFF CACHE BOOL "X11 not compatible with Fujitsu Sparc64 cross-compile, disabled." FORCE)

