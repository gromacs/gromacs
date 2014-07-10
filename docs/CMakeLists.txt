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

if(${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(MARKDOWN_CONFIGURE_IS_POSSIBLE off)
else()
    set(MARKDOWN_CONFIGURE_IS_POSSIBLE on)
endif()
find_package(Pandoc)

# Set up common infrastructure for configuring markdown at build time.
# Do replacement of CMake variables for version strings, etc. The use
# of configure-markdown.cmake defers until build time the
# configuration of markdown files, which could be faster for all the
# configurations that don't make the documentation even though it was
# possible, and helps avoid global re-configures if these files
# change.
set(SCRIPT_TO_CONFIGURE_MARKDOWN ${CMAKE_CURRENT_BINARY_DIR}/configure-markdown.cmake)
configure_file(configure-markdown.cmake.in
    ${SCRIPT_TO_CONFIGURE_MARKDOWN}
    @ONLY)

# Makes a custom command to configure a Markdown file found in the
# current source directory with the configure-markdown.make script
# produced above. The result is placed in the current binary directory
# for future use.
function(configure_markdown MARKDOWN_FILE)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${MARKDOWN_FILE}
        COMMAND ${CMAKE_COMMAND}
            -D FILE_TO_CONFIGURE=${CMAKE_CURRENT_SOURCE_DIR}/${MARKDOWN_FILE}
            -D CONFIGURED_FILE=${CMAKE_CURRENT_BINARY_DIR}/${MARKDOWN_FILE}
            -P ${SCRIPT_TO_CONFIGURE_MARKDOWN}
        DEPENDS
            ${SCRIPT_TO_CONFIGURE_MARKDOWN}
            ${CMAKE_CURRENT_SOURCE_DIR}/${MARKDOWN_FILE}
        COMMENT "Configuring Markdown"
        VERBATIM
        )
endfunction()

# Makes a custom command to make single-page HTML from Markdown. Takes
# a NAME argument for an output filename prefix, and a list of full
# paths to input files to concatenate with Pandoc into the HTML
# output.
function(make_markdown_html NAME)
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/${NAME}.html
        COMMAND
            ${PANDOC_EXECUTABLE} ${ARGN} -o ${HTML_OUTPUT_DIR}/${NAME}.html -s --toc --css buttondown.css
        DEPENDS ${ARGN}
        VERBATIM
        )
endfunction()

# Makes a custom command to make PDF from Markdown. Takes a NAME
# argument for an output filename prefix, and a list of full paths to
# input files to concatenate with Pandoc into the PDF.
function(make_markdown_pdf NAME)
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/${NAME}.pdf
        COMMAND
            ${PANDOC_EXECUTABLE} ${ARGN} -o ${HTML_OUTPUT_DIR}/${NAME}.pdf -s --toc
        DEPENDS ${ARGN}
        VERBATIM
        )
endfunction()

# function(make_markdown_multipage_html NAME)
#     # Make the multi-page HTML install guide
#
#     # TODO This is currently disabled, because the pandoc-specific
#     # buttondown.css doesn't work with the different kind of output
#     # makeinfo produces. When we understand better how we want to do
#     # generation, decide whether we want multi-page HTML output and
#     # how to make it work well.
#
#     add_custom_command(
#         OUTPUT ${HTML_OUTPUT_DIR}/${HTML_DIR}/index.html
#         COMMAND
#         ${PANDOC_EXECUTABLE} ${ARGN} -o ${NAME}.texi -s
#         COMMAND
#         ${MAKEINFO_EXECUTABLE} ${NAME}.texi --html -o ${HTML_OUTPUT_DIR}/${NAME} --css-ref buttondown.css
#         DEPENDS ${ARGN}
#         VERBATIM
#         )
# endfunction()

add_subdirectory(install-guide)
add_subdirectory(user-guide)
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

set(HTML_BUILD_IS_POSSIBLE OFF)
# We can only configure to build the webpage if the user asked for it,
# the build is outside of the source dir, and all the components can
# be built. There's no need to be talkative if we fail - most people
# never need to know.
if(GMX_BUILD_WEBPAGE AND
        GMX_BUILD_HELP AND
        NOT ${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR} AND
        MARKDOWN_CONFIGURE_IS_POSSIBLE AND
        MANUAL_BUILD_IS_POSSIBLE AND
        PANDOC_EXECUTABLE AND
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
    # the front page from index.md.
    if(BUILD_DOWNLOAD_SECTION)
        set(DOWNLOAD_SECTION
"# Downloads

*   Source code - [gromacs-${PROJECT_VERSION}.tar.gz](gromacs-${PROJECT_VERSION}.tar.gz)  
    (md5sum ${SOURCE_MD5SUM})  
    Other source code versions may be found at <ftp://ftp.gromacs.org/pub/gromacs/>

*   Regression tests - [regressiontests-${PROJECT_VERSION}.tar.gz](regressiontests-${PROJECT_VERSION}.tar.gz)  
    (md5sum ${REGRESSIONTESTS_MD5SUM})
")

        # Copy the source tarball to the webpage output
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E copy ${SOURCE_TARBALL} ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            VERBATIM
            )

        # Copy the regressiontests tarball to the webpage output
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            COMMAND ${CMAKE_COMMAND}
               -E copy ${REGRESSIONTESTS_TARBALL} ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            VERBATIM
            )
    else()
        set(DOWNLOAD_SECTION "")
    endif()

    # Put the CSS in the HTML output directory
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/buttondown.css
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/buttondown.css ${HTML_OUTPUT_DIR}/buttondown.css
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/buttondown.css
        VERBATIM
        )
    list(APPEND extra_webpage_dependencies ${HTML_OUTPUT_DIR}/buttondown.css)

    # Make the PDF reference guide
    # TODO Try to make the PDF arrive directly in ${HTML_OUTPUT_DIR}
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
        COMMAND ${CMAKE_COMMAND}
            -E remove -f ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
        COMMAND ${CMAKE_COMMAND}
            -E copy ${CMAKE_CURRENT_BINARY_DIR}/manual/gromacs.pdf ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
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

    # Add other dependencies for doing the webpage build from the real
    # tarball
    if(BUILD_DOWNLOAD_SECTION)
        list(APPEND extra_webpage_dependencies
            ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            )
    endif()

    configure_markdown(index.md)
    make_markdown_html(index ${CMAKE_CURRENT_BINARY_DIR}/index.md)

    # Add a top-level target for the others to hook onto
    add_custom_target(webpage
        DEPENDS
           ${HTML_OUTPUT_DIR}/index.html
           install-guide
           user-guide
           ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
           ${extra_webpage_dependencies}
        VERBATIM
        )
    add_dependencies(webpage webpage-html doc-all)
endif()
