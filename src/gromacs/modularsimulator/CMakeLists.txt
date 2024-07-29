#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
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

file(GLOB MODULARSIMULATOR_SOURCES *.cpp)

add_library(modularsimulator OBJECT ${MODULARSIMULATOR_SOURCES})

if (GMX_OPENMP)
    # Explicitly set properties for modular simulator module to compile with openmp
    set_target_properties(modularsimulator PROPERTIES COMPILE_OPTIONS $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_COMPILE_OPTIONS>)
endif ()

if (WIN32)
    gmx_target_warning_suppression(modularsimulator /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION)
    gmx_target_warning_suppression(modularsimulator /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
endif ()

target_include_directories(modularsimulator SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
gmx_target_compile_options(modularsimulator)
target_compile_definitions(modularsimulator PRIVATE HAVE_CONFIG_H)
# Should be possible to remove this when resolving #3290
target_include_directories(modularsimulator SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
target_compile_definitions(modularsimulator PRIVATE TMPI_USE_VISIBILITY)

# Source files have the following private external dependencies.
target_link_libraries(modularsimulator PRIVATE tng_io)

# Source files have the following private infrastructure dependencies.
target_link_libraries(modularsimulator PRIVATE common)

# Source files have the following private module dependencies.
# TODO: Explicitly link specific modules.
target_link_libraries(modularsimulator PRIVATE
                      legacy_api
                      legacy_modules
                      gmxlib
                      math
                      mdtypes
                      pbcutil
                      pulling
                      timing
                      tng_io
                      topology
                      utility
)

# Public interface for modules, including dependencies and interfaces
#target_include_directories(modularsimulator PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

if (BUILD_SHARED_LIBS)
    set_target_properties(modularsimulator PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
