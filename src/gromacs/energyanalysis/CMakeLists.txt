#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
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
add_library(energyanalysis OBJECT)
target_sources(energyanalysis PRIVATE
               energyterm.cpp
               )

if (BUILD_SHARED_LIBS)
    set_target_properties(energyanalysis PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

gmx_target_compile_options(energyanalysis)
target_compile_definitions(energyanalysis PRIVATE HAVE_CONFIG_H)

# Source files have the following dependencies on library infrastructure.
target_link_libraries(energyanalysis PRIVATE common)

# Public interface for modules, including dependencies and interfaces
#target_include_directories(energyanalysis PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(energyanalysis PUBLIC
                      fileio
                      legacy_api
                      utility
                      )

# Module dependencies
# energyanalysis interfaces convey transitive dependence on these modules.
target_link_libraries(energyanalysis PUBLIC
                      )
# Source files have the following private module dependencies.
#target_link_libraries(energyanalysis PRIVATE NOTHING)
# TODO: Explicitly link specific modules.
target_link_libraries(energyanalysis PRIVATE legacy_modules)

if(GMX_INSTALL_LEGACY_API)
    install(FILES
            DESTINATION include/gromacs/energyanalysis)
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
