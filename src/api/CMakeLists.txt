#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

# Activate targets for new C++ API components and docs.
if (GMXAPI)
    if (GMX_BUILD_MDRUN_ONLY)
        message(FATAL_ERROR "GMXAPI relies on libgromacs and is incompatible with GMX_BUILD_MDRUN_ONLY.")
    endif()
    if(NOT ${BUILD_SHARED_LIBS})
        # Note: this conditional should check for the existence of a libgromacs target supporting PIC
        # using the POSITION_INDEPENDENT_CODE property, but for now the only facility we have is the global
        # variable, BUILD_SHARED_LIBS.
        # TODO: gmxapi should gracefully build for dynamic linking or static linking for PIC or without.
        message(FATAL_ERROR "GMXAPI requires position-independent code. Set -DGMXAPI=OFF or -DBUILD_SHARED_LIBS=ON.")
    endif()
    add_subdirectory(cpp)
endif()

