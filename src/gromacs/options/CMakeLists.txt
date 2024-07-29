#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2010- The GROMACS Authors
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
add_library(options OBJECT
    abstractoption.cpp
    abstractsection.cpp
    basicoptions.cpp
    behaviorcollection.cpp
    filenameoption.cpp
    filenameoptionmanager.cpp
    options.cpp
    optionsassigner.cpp
    optionsection.cpp
    optionsvisitor.cpp
    timeunitmanager.cpp
    treesupport.cpp
    )

gmx_target_compile_options(options)
target_compile_definitions(options PRIVATE HAVE_CONFIG_H)

# Source files have the following private module dependencies.
target_link_libraries(options PRIVATE
                legacy_api
                legacy_modules
                utility
                )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(options PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(options PRIVATE common)

if (BUILD_SHARED_LIBS)  
    set_target_properties(options PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

