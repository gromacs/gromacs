#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2013- The GROMACS Authors
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

# Set up the module library
add_library(gmxana INTERFACE)

file(GLOB GMXANA_SOURCES *.cpp)

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${GMXANA_SOURCES} PARENT_SCOPE)

# Source files have the following private module dependencies.
target_link_libraries(gmxana PRIVATE
        )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(gmxana PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(gmxana PUBLIC
target_link_libraries(gmxana INTERFACE
        legacy_api
        )

# TODO: when gmxana is an OBJECT target
#target_link_libraries(gmxana PUBLIC legacy_api)
#target_link_libraries(gmxana PRIVATE common)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
