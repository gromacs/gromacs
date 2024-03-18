#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2018- The GROMACS Authors
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

set(TESTUTILS_DIR ${PROJECT_SOURCE_DIR}/src/testutils)
if (BUILD_TESTING)
    if(NOT GMX_DEVELOPER_BUILD)
        set(UNITTEST_TARGET_OPTIONS EXCLUDE_FROM_ALL)
    endif()
    include(${TESTUTILS_DIR}/TestMacros.cmake)
endif()

# The legacy public API headers may still be installed,
# but may be removed or changed without warning.
add_subdirectory(legacy)

# Activate targets for new C++ API components and docs.
if(GMX_NATIVE_WINDOWS)
    # GMXAPI has not been tested in Microsoft environments.
    # GMXAPI requires position-independent code
    set(_GMXAPI_DEFAULT OFF)
else()
    set(_GMXAPI_DEFAULT ${BUILD_SHARED_LIBS})
endif()
option(GMXAPI "Install GROMACS API." ${_GMXAPI_DEFAULT})
if (GMXAPI)
    if(NOT ${BUILD_SHARED_LIBS})
        # Note: this conditional should check for the existence of a libgromacs target supporting PIC
        # using the POSITION_INDEPENDENT_CODE property, but for now the only facility we have is the global
        # variable, BUILD_SHARED_LIBS.
        # TODO: gmxapi should gracefully build for dynamic linking or static linking for PIC or without.
        message(FATAL_ERROR "GMXAPI requires position-independent code. Set -DGMXAPI=OFF or -DBUILD_SHARED_LIBS=ON.")
    endif()
    add_subdirectory(gmxapi)
endif()

# Activate targets NBLIB
if(GMX_NATIVE_WINDOWS OR NOT BUILD_SHARED_LIBS)
    # NB-LIB has not been tested in Microsoft environments.
    # NB-LIB requires position-independent code
    set(_NBLIB_DEFAULT OFF)
else()
    set(_NBLIB_DEFAULT ${BUILD_SHARED_LIBS})
endif()
option(GMX_INSTALL_NBLIB_API "Install nblib headers" ${_NBLIB_DEFAULT})
if (GMX_INSTALL_NBLIB_API)
    if(NOT ${BUILD_SHARED_LIBS})
        message(FATAL_ERROR "NB-LIB requires position-independent code. Set -DGMX_INSTALL_NBLIB_API=OFF or -DBUILD_SHARED_LIBS=ON.")
    else()
        add_subdirectory(nblib)
    endif()
endif()

add_subdirectory(docs)
