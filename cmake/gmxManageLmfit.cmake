#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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

# Make a three-state enumeration, defaulting to 
gmx_option_multichoice(GMX_USE_LMFIT
    "How to handle the lmfit dependency of GROMACS"
    INTERNAL
    INTERNAL EXTERNAL NONE)
mark_as_advanced(GMX_USE_LMFIT)

# Make a fully functional lmfit library target that libgromacs can
# depend on regardless of how the user directed lmfit support and/or
# linking to work.
function(gmx_manage_lmfit)
    if(GMX_USE_LMFIT STREQUAL "INTERNAL")
        # Create an object library for the lmfit sources
        set(BUNDLED_LMFIT_DIR "${CMAKE_SOURCE_DIR}/src/external/lmfit")
        file(GLOB LMFIT_SOURCES ${BUNDLED_LMFIT_DIR}/*.cpp)
        add_library(lmfit_objlib OBJECT ${LMFIT_SOURCES})
        # Ensure that the objects can be used in both STATIC and SHARED
        # libraries.
        set_target_properties(lmfit_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

        # Create an INTERFACE (ie. fake) library for lmfit, that
        # libgromacs can depend on. The generator expression for the
        # target_sources expands to nothing when cmake builds the
        # export for libgromacs, so that it understands that we don't
        # install anything for this library - using plain source files
        # would not convey the right information.
        add_library(lmfit INTERFACE)
        target_sources(lmfit INTERFACE $<TARGET_OBJECTS:lmfit_objlib>)
        target_include_directories(lmfit SYSTEM INTERFACE $<BUILD_INTERFACE:${BUNDLED_LMFIT_DIR}>)
        # Add the lmfit interface library to the libgromacs Export name, even though
        # we will not be installing any content.
        install(TARGETS lmfit EXPORT libgromacs)

        set(HAVE_LMFIT 1 CACHE INTERNAL "Is lmfit found?")
    elseif(GMX_USE_LMFIT STREQUAL "EXTERNAL")
        # Find an external lmfit library.
        find_package(Lmfit ${GMX_LMFIT_REQUIRED_VERSION})
        if(NOT LMFIT_FOUND OR LMFIT_VERSION VERSION_LESS GMX_LMFIT_REQUIRED_VERSION)
            message(FATAL_ERROR "External lmfit >= ${GMX_LMFIT_REQUIRED_VERSION} could not be found, please adjust your pkg-config path to include the lmfit.pc file")
        endif()

        set(HAVE_LMFIT 1 CACHE INTERNAL "Is lmfit found?")
    else()
        # Create a dummy link target so the calling code doesn't need to know
        # whether lmfit support is being compiled.
        add_library(lmfit INTERFACE)
        # Add the lmfit interface library to the libgromacs Export name, even though
        # we will not be installing any content.
        install(TARGETS lmfit EXPORT libgromacs)

        set(HAVE_LMFIT 0 CACHE INTERNAL "Is lmfit found?")
    endif()
endfunction()
