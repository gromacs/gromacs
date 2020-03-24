#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020, by the GROMACS development team, led by
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
#
# \author Victor Holanda <victor.holanda@cscs.ch>
# \author Joe Jordan <ejjordan@kth.se>
# \author Prashanth Kanduri <kanduri@cscs.ch>
# \author Sebastian Keller <keller@cscs.ch>
#

file(GLOB NBLIB_SOURCES *.cpp)

add_library(nblib SHARED ${NBLIB_SOURCES})
gmx_target_compile_options(nblib)
target_compile_definitions(nblib PRIVATE HAVE_CONFIG_H)
target_include_directories(nblib SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
# Should be possible to change this when resolving #3290
target_include_directories(nblib SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)

###############################
# Install the public interface.
#
target_include_directories(nblib PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                           $<INSTALL_INTERFACE:include>
                           )
# If any item begins in a generator expression it must evaluate to a full path,
# so we cannot just use something like $<TARGET_PROPERTIES:nblibPublicHeaders,SOURCES>.
# Instead, we use a canonical list defined in the parent scope.
# install(DIRECTORY include/nblib
#         DESTINATION include)
# install(FILES
#     ${NBLIB_PUBLIC_HEADERS}
#     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nblib
#     COMPONENT libraries)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # Instruct a linking client to use its RPATH to resolve the nblib location.
    #
    # Explicitly specify library "install name" so that the correct loading
    # instruction is produced in client code. Client code should be able to find the
    # library relative to the client code RPATH. Without explicitly specifying,
    # INSTALL_NAME_DIR is inherited from the global CMAKE_INSTALL_NAME_DIR, which is
    # not appropriate for nblib if it uses an install name relative to the
    # executable_path or loader_path.
    set_target_properties(nblib PROPERTIES INSTALL_NAME_DIR "@rpath")
endif()

set_target_properties(nblib PROPERTIES
                      SOVERSION ${NBLIB_MAJOR}
                      VERSION ${NBLIB_RELEASE}
                      )

target_link_libraries(nblib PRIVATE libgromacs)


################################################
# Install and export nblib and Gromacs::nblib.
#
# Install the nblib target and simultaneously define the export target for
# which CMake will create a helper file. Specify the directory for clients to
# add to their include path to be able to `#include "nblib/some_header.h"`
install(TARGETS nblib
        EXPORT nblib
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )

# Create the CMake exports file to help other projects build against nblib
# as a CMake import target Gromacs::nblib.
install(EXPORT nblib
        NAMESPACE Gromacs::
        DESTINATION share/cmake/nblib/
        )
add_library(Gromacs::nblib ALIAS nblib )

include(CMakePackageConfigHelpers)

configure_package_config_file(
    cmake/nblib-config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/nblib-config.cmake"
    INSTALL_DESTINATION share/cmake/nblib/
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/nblib-config-version.cmake
    VERSION ${NBLIB_RELEASE}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/nblib-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/nblib-config.cmake
    DESTINATION share/cmake/nblib/
)

# We need a CMake target to provide the internal interface(s) of the nblib
# library implementation.
add_library(nblib-detail INTERFACE)
target_include_directories(nblib-detail
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

#if(BUILD_TESTING)
#    add_subdirectory(tests)
#    add_subdirectory(workflow/tests)
#endif()

#install(FILES
#        nbkerneloptions.h
#        DESTINATION include/nblib)

#set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${NBLIB_SOURCES} PARENT_SCOPE)

