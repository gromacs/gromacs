#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016,2018, by the GROMACS development team, led by
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

# Note that lmfit does not have a stable API, so GROMACS only supports
# the same version that it bundles.
set(GMX_LMFIT_REQUIRED_VERSION "7.0")

include(gmxOptionUtilities)

gmx_option_multichoice(GMX_USE_LMFIT
    "Use external lmfit instead of compiling the version bundled with GROMACS."
    INTERNAL
    INTERNAL EXTERNAL NONE)

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

function(manage_lmfit)
    if(GMX_USE_LMFIT STREQUAL "INTERNAL")
        set(BUNDLED_LMFIT_DIR "${CMAKE_SOURCE_DIR}/src/external/lmfit")
        file(GLOB LMFIT_SOURCES ${BUNDLED_LMFIT_DIR}/*.cpp)

        # Create a fake library for lmfit for libgromacs to depend on
        add_library(lmfit INTERFACE)
        target_sources(lmfit INTERFACE ${LMFIT_SOURCES})
        target_include_directories(lmfit INTERFACE
            $<BUILD_INTERFACE:${BUNDLED_LMFIT_DIR}>)
        target_link_libraries(libgromacs PRIVATE lmfit)

        set(HAVE_LMFIT_VALUE TRUE)
    elseif(GMX_USE_LMFIT STREQUAL "EXTERNAL")
        # Find an external lmfit library.
        find_package(Lmfit ${GMX_LMFIT_MINIMUM_REQUIRED_VERSION})
        if(NOT LMFIT_FOUND)
            message(FATAL_ERROR "External lmfit could not be found, please adjust your pkg-config path to include the lmfit.pc file")
        endif()
        target_link_libraries(libgromacs PRIVATE lmfit)
        set(HAVE_LMFIT_VALUE TRUE)
    else()
        set(HAVE_LMFIT_VALUE FALSE)
    endif()
    set(HAVE_LMFIT ${HAVE_LMFIT_VALUE} CACHE BOOL "Whether lmfit library support is enabled")
endfunction()
