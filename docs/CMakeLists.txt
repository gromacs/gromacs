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

# This directory provides a unified place for building all kinds of
# GROMACS documentation. This includes some "static" content (Doxygen
# code documentation, reference manual, install guide, old online HTML
# pages), and content generated from the gmx program for the various
# tools (man and HTML pages). It also provides the "webpage" target,
# that combines all of the above (except man pages in man format) into
# a form suitable for automated deployment to the GROMACS website.

add_subdirectory(man)
add_subdirectory(old-html)
add_subdirectory(doxygen)

# This would be prefereable to have as a sub-directory, but one
# add_custom_command() cannot depend on the output of another unless
# they are in the same file. Doing it like this is better than
# duplicating the logic to configure install-guide.md
include(BuildInstallGuide)

option(GMX_BUILD_WEBPAGE "Whether to try to configure to build the GROMACS static webpages" OFF)
mark_as_advanced(GMX_BUILD_WEBPAGE)

option(GMX_BUILD_MANUAL "Whether to try to configure to build the PDF manual" ${GMX_BUILD_WEBPAGE})
mark_as_advanced(GMX_BUILD_MANUAL)
if(GMX_BUILD_MANUAL)
    # Make sure we only do detection of manual-building dependencies
    # when the user opted in for that.
    add_subdirectory(manual)
endif()

find_package(Pandoc)

set(HTML_BUILD_IS_POSSIBLE OFF)
# We can only configure to build the webpage if the user asked for it,
# the build is outside of the source dir, and all the components can
# be built. There's no need to be talkative if we fail - most people
# never need to know.
if(GMX_BUILD_WEBPAGE AND
        GMX_BUILD_HELP AND
        NOT ${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR} AND
        INSTALL_GUIDE_BUILD_IS_POSSIBLE AND
        MANUAL_BUILD_IS_POSSIBLE AND
        PANDOC_EXECUTABLE AND
        DOXYGEN_EXECUTABLE AND
        DOXYGEN_MSCGEN_EXECUTABLE)
    set(HTML_BUILD_IS_POSSIBLE ON)
endif()

if(HTML_BUILD_IS_POSSIBLE)
    set(HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
    # Do replacement of CMake variables for version strings, etc.
    configure_file(configure-html.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/configure-html.cmake
        @ONLY)

    # For a real build of the webpage, the md5sum of the tarballs must
    # already be known, and so we may as well require that the real
    # build of the webpage take place from cmake run from the unpacked
    # tarball. Then, the *_MD5SUM and *_TARBALL variables will be able
    # to be set on the cmake command line (e.g. by a Jenkins job
    # configuration), and we can require that they are set. For local
    # building of the webpages (e.g. for debugging), the *_MD5SUM
    # variables need a fallback.
    if(NOT SOURCE_IS_SOURCE_DISTRIBUTION)
        foreach(VAR SOURCE_MD5SUM REGRESSIONTESTS_MD5SUM)
            if(NOT DEFINED ${VAR})
                set(${VAR} "unknown" CACHE STRING "Fallback value for output of md5sum")
                message(STATUS "Webpage build will not be fully featured, using fallback value for ${VAR}")
            endif()
        endforeach()
    else()
        foreach(VAR SOURCE_MD5SUM REGRESSIONTESTS_MD5SUM SOURCE_TARBALL REGRESSIONTESTS_TARBALL)
            if(NOT DEFINED ${VAR})
                message(FATAL_ERROR "The build of the webpage requires that ${VAR} is set in the cmake cache, e.g. on the command line")
            endif()
        endforeach()
    endif()

    # This defers until build time the configuration of
    # index.md, which could be faster
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/index.md
        COMMAND ${CMAKE_COMMAND}
            -P ${CMAKE_CURRENT_BINARY_DIR}/configure-html.cmake
        DEPENDS
            ${CMAKE_CURRENT_BINARY_DIR}/configure-html.cmake
            ${CMAKE_CURRENT_SOURCE_DIR}/index.md
        COMMENT "Configuring index.md"
        VERBATIM
        )

    # Put the CSS in the HTML output directory
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/buttondown.css
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/buttondown.css ${HTML_OUTPUT_DIR}/buttondown.css
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/buttondown.css
        VERBATIM
        )

    # Make the top-level index.html
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/index.html
        COMMAND ${PANDOC_EXECUTABLE} -o ${HTML_OUTPUT_DIR}/index.html index.md -s --css buttondown.css
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/index.md ${HTML_OUTPUT_DIR}/buttondown.css
        VERBATIM
        )

    # Make the HTML install guide
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/install-guide.html
        COMMAND
            ${PANDOC_EXECUTABLE} -t html -o ${HTML_OUTPUT_DIR}/install-guide.html install-guide.md -s --toc --css buttondown.css
        DEPENDS
            ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
            ${HTML_OUTPUT_DIR}/buttondown.css
        VERBATIM
        )

    # Make the PDF reference guide
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
        COMMAND ${CMAKE_COMMAND}
            -E copy ${CMAKE_CURRENT_BINARY_DIR}/manual/gromacs.pdf ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
        # UseLATEX.cmake makes a target called pdf, not ${CMAKE_CURRENT_BINARY_DIR}/manual/gromacs.pdf
        DEPENDS pdf
        VERBATIM
        )

    # Copy the source tarball to the webpage output, if available
    if(DEFINED SOURCE_TARBALL)
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E copy ${SOURCE_TARBALL} ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            VERBATIM
            )
    else()
        # Fallback, just touch the file
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E touch ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            VERBATIM
            )
    endif()

    # Copy the regressiontests tarball to the webpage output, if available
    if(DEFINED REGRESSIONTESTS_TARBALL)
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E copy ${REGRESSIONTESTS_TARBALL} ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            VERBATIM
            )
    else()
        # Fallback, just touch the file
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E touch ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            VERBATIM
            )
    endif()

    # TODO Move content from the "old" html output into the new user
    # guide, or delete, as appropriate.
    if(NOT SOURCE_IS_SOURCE_DISTRIBUTION)
        add_custom_target(webpage-html
            ${CMAKE_COMMAND} -E copy_directory old-html/final ${HTML_OUTPUT_DIR}
            )
        add_dependencies(webpage-html html)
    else()
        # In the source distribution, the html pages are already
        # built, so we can avoid making gmx via the html target
        add_custom_target(webpage-html
            ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/old-html/final ${HTML_OUTPUT_DIR}
            )
    endif()

    add_custom_target(webpage-doxygen
        ${CMAKE_COMMAND} -E make_directory ${HTML_OUTPUT_DIR}/doxygen
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/doxygen/index.html ${HTML_OUTPUT_DIR}/doxygen/
        COMMAND ${CMAKE_COMMAND} -E copy_directory doxygen/html-full ${HTML_OUTPUT_DIR}/doxygen/html-full
        COMMAND ${CMAKE_COMMAND} -E copy_directory doxygen/html-lib ${HTML_OUTPUT_DIR}/doxygen/html-lib
        COMMAND ${CMAKE_COMMAND} -E copy_directory doxygen/html-user ${HTML_OUTPUT_DIR}/doxygen/html-user
        )
    add_dependencies(webpage-doxygen doc-all)

    # Add a top-level target for the others to hook onto
    add_custom_target(webpage
        DEPENDS
           ${HTML_OUTPUT_DIR}/index.html
           ${HTML_OUTPUT_DIR}/install-guide.html
           ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
           ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
           ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
        VERBATIM
        )
    add_dependencies(webpage webpage-html webpage-doxygen)
endif()
