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

# This list file provides the Gromacs::gmxapi cmake target.

target_sources(gmxapi PRIVATE
               resourceassignment.cpp
               context.cpp
               exceptions.cpp
               gmxapi.cpp
               md.cpp
               mdmodule.cpp
               mdsignals.cpp
               session.cpp
               status.cpp
               system.cpp
               version.cpp
               workflow.cpp
               tpr.cpp
               )
set_target_properties(gmxapi PROPERTIES POSITION_INDEPENDENT_CODE ON)
gmx_target_compile_options(gmxapi)
target_compile_definitions(gmxapi PRIVATE HAVE_CONFIG_H)
target_include_directories(gmxapi SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
target_compile_definitions(gmxapi PRIVATE TMPI_USE_VISIBILITY)
# Should be possible to change this when resolving #3290
target_include_directories(gmxapi SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)

# Define implementation interface
target_include_directories(gmxapi PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           )

###############################
# Install the public interface.
#

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # Instruct a linking client to use its RPATH to resolve the libgmxapi location.
    #
    # Explicitly specify library "install name" so that the correct loading
    # instruction is produced in client code. Client code should be able to find the
    # library relative to the client code RPATH. Without explicitly specifying,
    # INSTALL_NAME_DIR is inherited from the global CMAKE_INSTALL_NAME_DIR, which is
    # not appropriate for libgmxapi if it uses an install name relative to the
    # executable_path or loader_path.
    set_target_properties(gmxapi PROPERTIES INSTALL_NAME_DIR "@rpath")
endif()

set_target_properties(gmxapi PROPERTIES
                      OUTPUT_NAME "gmxapi${GMX_LIBS_SUFFIX}"
                      SOVERSION ${GMXAPI_MAJOR}
                      VERSION ${GMXAPI_RELEASE}
                      )

target_link_libraries(gmxapi PRIVATE libgromacs)
target_link_libraries(gmxapi PRIVATE common)
# TODO: Remove `legacy_modules` once specific modules are explicitly linked.
target_link_libraries(gmxapi PRIVATE legacy_modules)
target_link_libraries(gmxapi PRIVATE math utility)


################################################
# Install and export gmxapi and Gromacs::gmxapi.
#
# Install the gmxapi target and simultaneously define the export target for
# which CMake will create a helper file. Specify the directory for clients to
# add to their include path to be able to `#include "gmxapi/some_header.h"`
install(TARGETS gmxapi
        EXPORT gmxapi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )

# Create the CMake exports file to help other projects build against libgmxapi
# as a CMake import target Gromacs::gmxapi.
install(EXPORT gmxapi
        NAMESPACE Gromacs::
        DESTINATION share/cmake/gmxapi/
        # If we use components for clients to ask for an MPI-aware version,
        # we might want to export the MPI-related details they need to help
        # with compatibility.
        # EXPORT_LINK_INTERFACE_LIBRARIES
        )
add_library(Gromacs::gmxapi ALIAS gmxapi )

include(CMakePackageConfigHelpers)

# WARNING: A hard-coded dependency on a specific gromacs configuration
# will likely be incompatible with a solution to #4334, but is a necessary
# workaround for #4563.
if(${CMAKE_VERSION} VERSION_EQUAL "3.24.0")
    set(
        _gmxapi_find_dependencies
        "set(gromacs${GMX_LIBS_SUFFIX}_ROOT \${PACKAGE_PREFIX_DIR})
find_dependency(gromacs${GMX_LIBS_SUFFIX} ${GMX_VERSION} CONFIG)"
    )
endif()
get_target_property(_mpi Gromacs::gmxapi MPI)
# TODO(#4334): Rename to gmxapiConfig.cmake to avoid globbing conflict when <name>.cmake includes <name>-*.cmake.
configure_package_config_file(
    cmake/gmxapi-config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config.cmake"
    INSTALL_DESTINATION share/cmake/gmxapi/
)
unset(_mpi)
unset(_gmxapi_find_dependencies)
get_target_property(GMXAPI_RELEASE gmxapi RELEASE)
# TODO(#4334): Rename to gmxapiConfigVersion.cmake to avoid globbing conflict when <name>.cmake includes <name>-*.cmake.
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config-version.cmake
    VERSION ${GMXAPI_RELEASE}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config.cmake
    DESTINATION share/cmake/gmxapi/
)

# We need a CMake target to provide the internal interface(s) of the gmxapi
# library implementation.
add_library(gmxapi-detail INTERFACE)
target_include_directories(gmxapi-detail
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

if (BUILD_TESTING AND GMX_BUILD_UNITTESTS)
    add_subdirectory(tests)
    add_subdirectory(workflow/tests)
endif()
