#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018, by the GROMACS development team, led by
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

# This list file provides the gmxapi::gmxapi cmake module. gmxapi::tests is the lone
# submodule, but subdirectories include the public headers.

##########################
# Set up public interface.
#
# A product in the parent (src/api/) CMake scope is a variable named
# GMXAPI_PUBLIC_BUILD_HEADERS containing the full paths to the public
# headers that are being installed.

add_subdirectory(include)

# The include directory should be mostly empty so that we can use it internally as
# the public interface include directory during build and testing.
configure_file(include/version.h.in include/gmxapi/version.h)

# Explicitly declare the public interface. This variable has local scope. Outside
# of this file, the SOURCES property of the gmxapiPublicHeaders target can be
# inspected or, for full paths in the build environment, the GMXAPI_PUBLIC_BUILD_HEADERS
# variable. See below.
set(_gmxapi_interface_headers
    include/gmxapi/context.h
    include/gmxapi/exceptions.h
    include/gmxapi/gmxapi.h
    include/gmxapi/gromacsfwd.h
    include/gmxapi/md.h
    include/gmxapi/session.h
    include/gmxapi/status.h
    include/gmxapi/system.h
    include/gmxapi/version.h)

################################################################################
# Create a helpful GMXAPI_PUBLIC_BUILD_HEADERS variable in the src/api/ CMake scope.
#
# Collect the public sources. Let CMake find the full path and we can use the
# SOURCES property and LOCATION properties to get the full file paths for
# installation and for other targets, such as the documentation builder inputs.
# Otherwise, this object target is not used by the build system.
#
# Alternative: write a list of the public headers to a file with a custom command.
# I think there is a CMake-ism to find the output of a custom command from a
# different directory.
# Custom commands can use generator expressions. Then the list of files can be
# loaded into a variable with file(STRINGS filename) and maybe iterative
# file(TO_NATIVE_PATH ....
# file(GENERATRE ...) can also process generator expressions, but the output
# file will not be available at generation time, and so will not be available
# to other file() commands.
# I can't think of any other way
# to use the $<TARGET_PROPERTY:gmxapi_public,SOURCES> to generate the
# list of officially public headers without duplication.
#
# This automation may be more confusing than helpful. We could just explicitly
# list the static and generated headers separately and convert them to the
# helper variable with string processing.

add_library(gmxapiPublicHeaders OBJECT ${_gmxapi_interface_headers})
target_include_directories(gmxapiPublicHeaders PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/include
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/include
                           )

# Create a variable in src/api/ scope with the full paths to the public interface.
# Reset the relative paths in _gmxapi_public_build_headers and get absolute paths.
set(_gmxapi_public_build_headers "")
get_target_property(_names gmxapiPublicHeaders SOURCES)
foreach(_name ${_names})
    get_source_file_property(_path ${_name} LOCATION)
    list(APPEND _gmxapi_public_build_headers ${_path})
endforeach()
# Note: GMXAPI_PUBLIC_BUILD_HEADERS is _only_ available in the parent scope. In
# the current scope, it is empty.
set(GMXAPI_PUBLIC_BUILD_HEADERS ${_gmxapi_public_build_headers} PARENT_SCOPE)
unset(_name)
unset(_names)

###############################
# Install the public interface.
#
# If any item begins in a generator expression it must evaluate to a full path,
# so we can't just use $<TARGET_PROPERTIES:gmxapiPublicHeaders,SOURCES>.
# Also, GMXAPI_PUBLIC_BUILD_HEADERS is not in this scope, so use _gmxapi_public_build_headers
install(FILES ${_gmxapi_public_build_headers}
        DESTINATION include/gmxapi)

unset(_sources_list)

##############################
# Set up gmxapi library target
#
# Define the gmxapi library target. Dependent targets can use properties of the target
# for include directories and header files.
# We explicitly provide the public headers to help development tools understand
# the environment (include paths, etc.) in which the headers will be used. CMake
# is smart enough not to set the HEADER_FILE_ONLY property and not compile *.h
# We do not link against gmxapiPublicHeaders because it adds unnecessary extra
# indirection and we want to be clear about all where the properties of gmxapi
# come from.

add_library(gmxapi SHARED ${_gmxapi_interface_headers})

# Define public interface. Make sure targets linking against `gmxapi` in the build
# system don't accidentally have the implementation headers (this directory))
# in a default include path.
target_include_directories(gmxapi PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/include
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/include
    $<INSTALL_INTERFACE:include>
    )

# Define private interface in build environment. (Not installed.)
# Unit tests for implementation details can use this target.
add_library(gmxapiImpl INTERFACE)
target_include_directories(gmxapiImpl INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           )
add_library(Gromacs::gmxapiImpl ALIAS gmxapiImpl)

# Define implementation interface
target_include_directories(gmxapi PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}>/include
                           ${CMAKE_CURRENT_BINARY_DIR}>/include
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           )

# libgmxapi source files.
# Listing all source files explicitly helps some code introspection tools
# and IDEs work better.

# Private / implementation headers in src/api/cpp
target_sources(gmxapi PRIVATE
               md-impl.h
               mdsignals-impl.h
               session-impl.h
               sessionresources-impl.h
               system-impl.h
               workflow-impl.h)

target_sources(
    gmxapi PRIVATE
    context.cpp
    gmxapi.cpp
    md.cpp
    mdmodule.cpp
    mdsignals.cpp
    outputstream.cpp
    session.cpp
    status.cpp
    system.cpp
    version.cpp
    workflow.cpp
)

# Help development tools to understand how include paths will be resolved.
target_sources(gmxapi PRIVATE include/version.h.in)
# This property should be set automatically, but this line is left here as a
# reminder in case version.in.h gets renamed to version.h.in at some point.
set_source_files_properties(include/version.h.in PROPERTIES HEADER_FILE_ONLY ON)

# Ref. https://cmake.org/Wiki/CMake_RPATH_handling
# use, i.e. don't skip the full RPATH for the build tree
set_target_properties(gmxapi PROPERTIES SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set_target_properties(gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)

# May duplicate link path when installed in LIB_INSTALL_DIR with libgromacs
# but without it, targets referencing libgromacs symbols when importing only this target will
# not find libgromacs the way OS X has set the RPATH to @executable_path=../lib
# To set the @rpath to resolve to the directory in which the libgmxapi shared
# object is living, we would set it to @loader_path, but I don't know of a CMake
# directive that produces that effect. We may just need a platform check.
# I suppose the question is what/how/why does CMake set rpath to executable_path instead
# of rpath...
set_target_properties(gmxapi PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

set_target_properties(gmxapi PROPERTIES
                      INSTALL_RPATH_USE_LINK_PATH TRUE
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
    PATH_VARS CMAKE_INSTALL_LIBDIR
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

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
