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
add_library(fileio INTERFACE)
file(GLOB FILEIO_SOURCES *.cpp)

if(GMX_USE_PLUGINS)
    add_library(vmddlopen OBJECT ${CMAKE_SOURCE_DIR}/src/external/vmd_molfile/vmddlopen.cpp)
    gmx_target_compile_options(vmddlopen)
    if (WIN32)
        gmx_target_warning_suppression(vmddlopen /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
    endif()
    list(APPEND libgromacs_object_library_dependencies vmddlopen)
    set(libgromacs_object_library_dependencies ${libgromacs_object_library_dependencies} PARENT_SCOPE)
else()
    # Remove vmdio.cpp from sources since we do not have plugin support
    list(FILTER FILEIO_SOURCES EXCLUDE REGEX ".*vmdio.cpp$")
endif()

# Source files have the following private module dependencies.
target_link_libraries(fileio PRIVATE
#                      gmxlib
#                      math
#                      mdtypes
#                      tng_io
                      )

# Only needed for XDR routines that cast xdr_char and the like to xdrproc_t.
# TODO: Convert to gmx_target_warning_suppression when fileio in an OBJECT target
gmx_target_interface_warning_suppression(fileio "-Wno-cast-function-type-strict" HAS_WARNING_NO_CAST_FUNCTION_TYPE_STRICT)

# Public interface for modules, including dependencies and interfaces
#target_include_directories(fileio PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(fileio PUBLIC
target_link_libraries(fileio INTERFACE
                      legacy_api
                      math
                      )

# TODO: when fileio is an OBJECT target
#target_link_libraries(fileio PUBLIC legacy_api)
#target_link_libraries(fileio PRIVATE common)

# Module dependencies
# fileio interfaces convey transitive dependence on these modules.
#target_link_libraries(fileio PUBLIC
target_link_libraries(fileio INTERFACE
                      utility
                      )
# Source files have the following private module dependencies.
#target_link_libraries(fileio PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(fileio PRIVATE legacy_modules)

target_sources(libgromacs PRIVATE ${FILEIO_SOURCES})

if (BUILD_TESTING)
     add_subdirectory(tests)
endif()
