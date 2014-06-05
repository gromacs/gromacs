#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

set(INSTALL_GUIDE_BUILD_IS_POSSIBLE OFF)
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    # in-source build would not be stable
    message(STATUS "Cannot build the install guide when building in the source directory")
else()
    find_package(Pandoc)
    if(PANDOC_EXECUTABLE)
        set(INSTALL_GUIDE_BUILD_IS_POSSIBLE ON)
    endif()
endif()

if(INSTALL_GUIDE_BUILD_IS_POSSIBLE)
    # Do replacement of CMake variables for version strings, etc.
    configure_file(configure-install-guide.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/configure-install-guide.cmake
        @ONLY)

    # This defers until build time the configuration of
    # install-guide.md, which could be faster
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
        COMMAND ${CMAKE_COMMAND}
            -P ${CMAKE_CURRENT_BINARY_DIR}/configure-install-guide.cmake
        DEPENDS
            ${CMAKE_CURRENT_BINARY_DIR}/configure-install-guide.cmake
            ${CMAKE_CURRENT_SOURCE_DIR}/install-guide.md
        COMMENT "Configuring install guide"
        VERBATIM
        )

    # Make the HTML install guide
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/install-guide.html
        COMMAND pandoc -t html -o ${CMAKE_CURRENT_BINARY_DIR}/install-guide.html install-guide.md -s --toc
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
        VERBATIM
        )

    # Make the INSTALL file for CPack for the tarball
    add_custom_command(
        OUTPUT ${CMAKE_BINARY_DIR}/INSTALL
        COMMAND pandoc -t plain -o ../INSTALL install-guide.md
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
        VERBATIM
        )

    # Add a top-level target for the others to hook onto
    add_custom_target(install-guide
        DEPENDS
           ${CMAKE_CURRENT_BINARY_DIR}/install-guide.html
           ${CMAKE_BINARY_DIR}/INSTALL
        VERBATIM
        )
endif()
