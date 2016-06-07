#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016, by the GROMACS development team, led by
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

set(GMX_LMFIT_MINIMUM_REQUIRED_VERSION "6.1")
set(GMX_BUNDLED_LMFIT_DIR "${CMAKE_SOURCE_DIR}/src/external/lmfit")

option(GMX_EXTERNAL_LMFIT "Use external lmfit instead of compiling the version bundled with GROMACS." OFF)
mark_as_advanced(GMX_EXTERNAL_LMFIT)

macro(manage_lmfit)
    if(GMX_EXTERNAL_LMFIT)
        # Find an external lmfit library.
        find_package(Lmfit ${GMX_LMFIT_MINIMUM_REQUIRED_VERSION})
        if(NOT LMFIT_FOUND)
            message(FATAL_ERROR "External lmfit could not be found, please adjust your pkg-config path to include the lmfit.pc file")
        endif()
    endif()
endmacro()

macro(get_lmfit_properties LMFIT_SOURCES_VAR LMFIT_LIBRARIES_VAR LMFIT_INCLUDE_DIR_VAR LMFIT_INCLUDE_DIR_ORDER_VAR)
    if (GMX_EXTERNAL_LMFIT)
        set(${LMFIT_INCLUDE_DIR_VAR} ${LMFIT_INCLUDE_DIR})
        set(${LMFIT_INCLUDE_DIR_ORDER_VAR} "AFTER")
        set(${LMFIT_SOURCES_VAR} "")
        set(${LMFIT_LIBRARIES_VAR} ${LMFIT_LIBRARIES})
    else()
        set(${LMFIT_INCLUDE_DIR_VAR} ${GMX_BUNDLED_LMFIT_DIR})
        set(${LMFIT_INCLUDE_DIR_ORDER_VAR} "BEFORE")
        file(GLOB ${LMFIT_SOURCES_VAR} ${GMX_BUNDLED_LMFIT_DIR}/*.cpp)
        set(${LMFIT_LIBRARIES_VAR} "")
    endif()
endmacro()

manage_lmfit()

