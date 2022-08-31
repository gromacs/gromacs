#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2017- The GROMACS Authors
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
add_library(taskassignment OBJECT)
target_sources(taskassignment PRIVATE
               decidegpuusage.cpp
               decidesimulationworkload.cpp
               findallgputasks.cpp
               reportgpuusage.cpp
               resourcedivision.cpp
               taskassignment.cpp
               usergpuids.cpp)
if (BUILD_SHARED_LIBS)
    set_target_properties(taskassignment PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Note that this is a higher-level module that should not need
# special compilation to suit e.g. GPU or MPI dependencies.
# If you find you want to do that, consider a preliminary
# refactoring.

# Source files have the following dependencies on library infrastructure.
target_link_libraries(taskassignment PRIVATE common)

# Public interface for modules, including dependencies and interfaces
target_include_directories(taskassignment PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(taskassignment PRIVATE
                      legacy_api
                      legacy_modules
                      math
                      pulling
                      timing
                      topology
                      utility
                      )

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
