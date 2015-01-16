#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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
# a form suitable for automated deployment to the GROMACS website. It
# also provides the INSTALL file for the tarball.
#
# All of the markdown content is configured, and we'd like to do that
# at build time rather than configure time (for speed, when not
# building markdown content). Also, the way they should be configured
# varies with whether the source is a tarball or repo, and which file
# is being configured. So several *_IS_POSSIBLE variables are used to
# direct the configure-time logic so that all appropriate variables
# are set by the time the configure-markdown.cmake.in file is
# configured, so that later it can do the configuration of all the
# markdown files and the right thing will happen in each case.

# Even if we aren't going to make the full webpage, set up to put all
# the documentation output in the same place, for convenience
set(HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
file(MAKE_DIRECTORY ${HTML_OUTPUT_DIR})

find_package(Sphinx 1.2.3 COMPONENTS pygments)

# Set up common infrastructure for configuring sphinx files at build
# time.  Do replacement of CMake variables for version strings,
# etc. The use of configure-sphinx.cmake defers until build time the
# configuration of sphinx config and .rst files, which could be faster
# for all the configurations that don't make the documentation even
# though it was possible, and helps avoid global re-configures if
# these files change.
set(SCRIPT_TO_CONFIGURE_SPHINX ${CMAKE_CURRENT_BINARY_DIR}/configure-sphinx.cmake)
configure_file(configure-sphinx.cmake.in
    ${SCRIPT_TO_CONFIGURE_SPHINX}
    @ONLY)

# Makes a custom command to configure a Sphinx file found in the
# current source directory with the configure-sphinx.cmake script
# produced above. The result is placed in the current binary directory
# for future use.
set(SPHINX_INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/sphinx-input)
function(configure_sphinx SPHINX_FILE)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${SPHINX_FILE}
        COMMAND ${CMAKE_COMMAND}
            -D FILE_TO_CONFIGURE=${CMAKE_CURRENT_SOURCE_DIR}/${SPHINX_FILE}
            -D CONFIGURED_FILE=${SPHINX_INPUT_DIR}/${SPHINX_FILE}
            -P ${SCRIPT_TO_CONFIGURE_SPHINX}
        DEPENDS
            ${SCRIPT_TO_CONFIGURE_SPHINX}
            ${CMAKE_CURRENT_SOURCE_DIR}/${SPHINX_FILE}
        COMMENT "Configuring Sphinx file ${SPHINX_FILE}"
        VERBATIM
        )
endfunction()

add_subdirectory(man)
add_subdirectory(old-html)
add_subdirectory(doxygen)

option(GMX_BUILD_WEBPAGE "Whether to try to configure to build the GROMACS static webpages" OFF)
mark_as_advanced(GMX_BUILD_WEBPAGE)

option(GMX_BUILD_MANUAL "Whether to try to configure to build the PDF manual" ${GMX_BUILD_WEBPAGE})
mark_as_advanced(GMX_BUILD_MANUAL)
if(GMX_BUILD_MANUAL)
    # Make sure we only do detection of manual-building dependencies
    # when the user opted in for that.
    add_subdirectory(manual)
endif()

if(NOT BUILD_IS_INSOURCE AND
        SPHINX_FOUND)
    # Configure files that are inputs to Sphinx, either for the
    # webpage build or (in a few cases) for building the INSTALL file
    # for the tarball.
    list(APPEND RST_FILES
        "index.rst"
        "download.rst"
        "online.rst"
        "install-guide/main.rst"
        "user-guide/main.rst"
        "user-guide/introduction.rst"
        "user-guide/getting-started.rst"
        "user-guide/file-formats.rst"
        "user-guide/tools.rst"
        "conf.py")
    foreach(file ${RST_FILES})
        configure_sphinx(${file})
        list(APPEND BUILD_DIR_RST_FILES "${CMAKE_CURRENT_BINARY_DIR}/${file}")
    endforeach()
    set(SPHINX_FILES_WERE_CONFIGURED 1)
    # Make the INSTALL file for CPack for the tarball. This gets put
    # into the tarball via the CPack rules in the top-level
    # CMakeLists.txt. Note that because sphinx-build relies on
    # recursive globbing of .rst files, TARBALL_INSTALL_GUIDE_DIR must
    # be distinct from SPHINX_INPUT_DIR.
    set(TARBALL_INSTALL_GUIDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/install-guide-for-tarball")
    set(TARBALL_INSTALL_GUIDE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/install-guide/final")
    add_custom_command(
        OUTPUT final/INSTALL
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARBALL_INSTALL_GUIDE_DIR}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARBALL_INSTALL_GUIDE_OUTPUT_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy ${SPHINX_INPUT_DIR}/install-guide/main.rst ${TARBALL_INSTALL_GUIDE_DIR}/index.rst
        COMMAND ${CMAKE_COMMAND} -E copy ${SPHINX_INPUT_DIR}/conf.py ${TARBALL_INSTALL_GUIDE_DIR}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TARBALL_INSTALL_GUIDE_OUTPUT_DIR}
        COMMAND ${SPHINX_EXECUTABLE} -q -b text ${TARBALL_INSTALL_GUIDE_DIR} ${TARBALL_INSTALL_GUIDE_OUTPUT_DIR}
        COMMAND ${CMAKE_COMMAND} -E rename ${TARBALL_INSTALL_GUIDE_OUTPUT_DIR}/index.txt ${TARBALL_INSTALL_GUIDE_OUTPUT_DIR}/INSTALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/install-guide/main.rst ${CMAKE_CURRENT_BINARY_DIR}/conf.py
        VERBATIM
        )

    # Add a convenience top-level target for the Jenkins tarball build
    add_custom_target(install-guide
        DEPENDS final/INSTALL
        VERBATIM
        )
    gmx_cpack_add_generated_source_directory(final DESTINATION /)
endif()

set(HTML_BUILD_IS_POSSIBLE OFF)
# We can only configure to build the webpage if the user asked for it,
# the build is outside of the source dir, and all the components can
# be built. There's no need to be talkative if we fail - most people
# never need to know.
if(GMX_BUILD_WEBPAGE AND
        GMX_BUILD_HELP AND
        SPHINX_FILES_WERE_CONFIGURED AND
        DOXYGEN_EXECUTABLE AND
        DOXYGEN_MSCGEN_EXECUTABLE)
    set(HTML_BUILD_IS_POSSIBLE ON)
endif()

if(HTML_BUILD_IS_POSSIBLE)
    # For a real build of the webpage, the md5sum of the tarballs must
    # already be known, and so we may as well require that the real
    # build of the webpage take place from cmake run from the unpacked
    # tarball. Then, the *_MD5SUM and *_TARBALL variables will be able
    # to be set on the cmake command line (e.g. by a Jenkins job
    # configuration), and we can require that they are set. For local
    # building of the webpages (e.g. for debugging), those variables
    # can be left unset, and if so, the download section will not be
    # constructed.
    if(NOT SOURCE_IS_SOURCE_DISTRIBUTION)
        if (SOURCE_TARBALL AND SOURCE_MD5SUM AND
                REGRESSIONTESTS_TARBALL AND REGRESSIONTESTS_MD5SUM)
            set(BUILD_DOWNLOAD_SECTION on)
        else()
            set(BUILD_DOWNLOAD_SECTION off)
        endif()
    else()
        foreach(VAR SOURCE_MD5SUM REGRESSIONTESTS_MD5SUM SOURCE_TARBALL REGRESSIONTESTS_TARBALL)
            if(NOT DEFINED ${VAR})
                message(FATAL_ERROR "The build of the webpage requires that ${VAR} is set in the cmake cache, e.g. on the CMake command line")
            endif()
        endforeach()
        set(BUILD_DOWNLOAD_SECTION on)
    endif()

    # If building the webpage from the repo, then tarballs may not
    # exist, and if so, it would not make sense to build that part of
    # the front page from index.rst.
    if(BUILD_DOWNLOAD_SECTION)
        # Copy the source tarball to the webpage output
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/gromacs-${GMX_VERSION_STRING}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E copy ${SOURCE_TARBALL} ${HTML_OUTPUT_DIR}/gromacs-${GMX_VERSION_STRING}.tar.gz
            VERBATIM
            )

        # Copy the regressiontests tarball to the webpage output
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/regressiontests-${GMX_VERSION_STRING}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E copy ${REGRESSIONTESTS_TARBALL} ${HTML_OUTPUT_DIR}/regressiontests-${GMX_VERSION_STRING}.tar.gz
            VERBATIM
            )

        # Add other dependencies for doing the webpage build from the real
        # tarball
        list(APPEND extra_webpage_dependencies
            ${HTML_OUTPUT_DIR}/gromacs-${GMX_VERSION_STRING}.tar.gz
            ${HTML_OUTPUT_DIR}/regressiontests-${GMX_VERSION_STRING}.tar.gz
            )

        # Get Sphinx to do the right thing in index.rst
        set(BUILD_DOWNLOAD_SECTION_TAG "-t is_webpage")
    endif()

    # Make the PDF reference guide
    # TODO Try to make the PDF arrive directly in ${HTML_OUTPUT_DIR}
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/manual-${GMX_VERSION_STRING}.pdf
        COMMAND ${CMAKE_COMMAND}
            -E remove -f ${HTML_OUTPUT_DIR}/manual-${GMX_VERSION_STRING}.pdf
        COMMAND ${CMAKE_COMMAND}
            -E copy ${CMAKE_CURRENT_BINARY_DIR}/manual/gromacs.pdf ${HTML_OUTPUT_DIR}/manual-${GMX_VERSION_STRING}.pdf
        # UseLATEX.cmake makes a target called pdf, not ${CMAKE_CURRENT_BINARY_DIR}/manual/gromacs.pdf
        DEPENDS pdf
        VERBATIM
        )

    # TODO Move content from the "old" html output into the new user
    # guide, or delete, as appropriate.
    if(NOT SOURCE_IS_SOURCE_DISTRIBUTION)
        # TODO If content remains here once the user guide is in
        # decent shape, try to make the generated HTML arrive directly
        # in ${HTML_OUTPUT_DIR}
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

    # The Doxygen configuration in doxygen/Doxyfile-common.cmakein
    # makes all the Doxygen output directly in
    # ${HTML_OUTPUT_DIR}/doxygen (and makes the directory if it needs
    # to).

    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

    add_custom_target(webpage-sphinx ALL
        COMMAND
            cmake -E make_directory ${SPHINX_INPUT_DIR}/_static
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b html
            -d "${SPHINX_CACHE_DIR}"
            "${BUILD_DOWNLOAD_SECTION_TAG}"
            "${SPHINX_INPUT_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}/html"
        DEPENDS
            ${BUILD_DIR_RST_FILES}
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building HTML documentation with Sphinx"
        )

    # Add a top-level target for the others to hook onto
    add_custom_target(webpage
        DEPENDS
           ${HTML_OUTPUT_DIR}/index.html
           webpage-sphinx
           ${HTML_OUTPUT_DIR}/manual-${GMX_VERSION_STRING}.pdf
           ${extra_webpage_dependencies}
        VERBATIM
        )
    add_dependencies(webpage webpage-html doxygen-all)
endif()
