#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019, by the GROMACS development team, led by
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

# This list file provides the Gromacs::gmxapi cmake module.

##########################
# Set up public interface.
#
# The parent (src/api/) CMake scope provides a variable named
# GMXAPI_PUBLIC_HEADERS containing the full paths to the public
# headers that are being installed. The headers are segregated into a
# subdirectory here so that their build-time include directory path
# does not expose lower level headers.

add_subdirectory(include)

# The include directory should be mostly empty so that we can use it internally as
# the public interface include directory during build and testing.
configure_file(include/version.h.in include/gmxapi/version.h)

add_library(gmxapi SHARED
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
            workflow.cpp)

# Define public interface. Make sure targets linking against `gmxapi` in the build
# system don't accidentally have the implementation headers (this directory))
# in a default include path.
target_include_directories(gmxapi PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                           $<INSTALL_INTERFACE:include>
                           )
# Define implementation interface
target_include_directories(gmxapi PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           )

###############################
# Install the public interface.
#
# If any item begins in a generator expression it must evaluate to a full path,
# so we can't just use something like $<TARGET_PROPERTIES:gmxapiPublicHeaders,SOURCES>.
# Instead, we use a canonical list defined in the parent scope.
install(DIRECTORY include/gmxapi
        DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/gmxapi/version.h
        DESTINATION include/gmxapi)

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
                      SOVERSION ${GMXAPI_MAJOR}
                      VERSION ${GMXAPI_RELEASE}
                      )

target_link_libraries(gmxapi PRIVATE libgromacs)


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
        )
add_library(Gromacs::gmxapi ALIAS gmxapi )

include(CMakePackageConfigHelpers)

configure_package_config_file(
    cmake/gmxapi-config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config.cmake"
    INSTALL_DESTINATION share/cmake/gmxapi/
)
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

if(BUILD_TESTING)
    add_subdirectory(tests)
    add_subdirectory(workflow/tests)
endif()
