#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014 by the GROMACS development team.
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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

add_library(analysisdata INTERFACE)
file(GLOB ANALYSISDATA_SOURCES *.cpp modules/*.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${ANALYSISDATA_SOURCES} PARENT_SCOPE)

if(GMX_INSTALL_LEGACY_API)
    install(FILES
            abstractdata.h
            analysisdata.h
            arraydata.h
            dataframe.h
            datamodule.h
            DESTINATION include/gromacs/analysisdata)
endif()

# Source files have the following private module dependencies.
target_link_libraries(analysisdata PRIVATE
#                      gmxlib
#                      math
#                      mdtypes
#                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(analysisdata PUBLIC
target_include_directories(analysisdata INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(analysisdata PUBLIC
target_link_libraries(analysisdata INTERFACE
                      legacy_api
                      )

# TODO: when analysisdata is an OBJECT target
#target_link_libraries(analysisdata PUBLIC legacy_api)
#target_link_libraries(analysisdata PRIVATE common)

# Module dependencies
# analysisdata interfaces convey transitive dependence on these modules.
#target_link_libraries(analysisdata PUBLIC
target_link_libraries(analysisdata INTERFACE
                      utility
                      )
# Source files have the following private module dependencies.
#target_link_libraries(analysisdata PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(analysisdata PRIVATE legacy_modules)

add_subdirectory(modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
