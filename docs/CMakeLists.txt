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

add_subdirectory(man)
add_subdirectory(old-html)
add_subdirectory(doxygen)

if(${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(MARKDOWN_CONFIGURE_IS_POSSIBLE off)
else()
    set(MARKDOWN_CONFIGURE_IS_POSSIBLE on)
endif()

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
endif()

# Note that the install-guide target, etc. can still be built
# independent of the webpage target.
if(MARKDOWN_CONFIGURE_IS_POSSIBLE AND PANDOC_EXECUTABLE)
    # Do replacement of CMake variables for version strings, etc.  The
    # use of configure-markdown.cmake defers until build time the
    # configuration of markdown files, which could be faster for all
    # the configurations that don't make the documentation even though
    # it was possible, and helps avoid global re-configures if these
    # files change.
    configure_file(configure-markdown.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/configure-markdown.cmake
        @ONLY)

    # Configure the markdown at build time
    foreach(filename_prefix index install-guide user-guide)
        add_custom_command(
            OUTPUT ${filename_prefix}.md
            COMMAND ${CMAKE_COMMAND} -DFILENAME_PREFIX_TO_CONFIGURE=${filename_prefix}
                -P ${CMAKE_CURRENT_BINARY_DIR}/configure-markdown.cmake
            DEPENDS
                ${CMAKE_CURRENT_BINARY_DIR}/configure-markdown.cmake
                ${CMAKE_CURRENT_SOURCE_DIR}/${filename_prefix}.md
            COMMENT "Configuring Markdown"
            VERBATIM
        )
    endforeach()

    # Put the CSS in the HTML output directory
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/buttondown.css
        COMMAND
            ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/buttondown.css ${HTML_OUTPUT_DIR}/buttondown.css
        DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/buttondown.css
        VERBATIM
        )

    # Make the top-level index.html
    add_custom_command(
        OUTPUT ${HTML_OUTPUT_DIR}/index.html
        COMMAND
            ${PANDOC_EXECUTABLE} -o ${HTML_OUTPUT_DIR}/index.html index.md -s --css buttondown.css
        DEPENDS
            ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
        VERBATIM
        )

    # Make the single-page HTML and PDF versions of install and user guides
    foreach(type install user)
        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/${type}-guide.html
            COMMAND
                ${PANDOC_EXECUTABLE} ${type}-guide.md -o ${HTML_OUTPUT_DIR}/${type}-guide.html -s --toc --css buttondown.css
            DEPENDS
                ${CMAKE_CURRENT_BINARY_DIR}/${type}-guide.md
            VERBATIM
            )

        add_custom_command(
            OUTPUT ${HTML_OUTPUT_DIR}/${type}-guide.pdf
            COMMAND
                ${PANDOC_EXECUTABLE} ${type}-guide.md -o ${HTML_OUTPUT_DIR}/${type}-guide.pdf -s --toc
            DEPENDS
                ${CMAKE_CURRENT_BINARY_DIR}/${type}-guide.md
            VERBATIM
            )
    endforeach()

    # Make the multi-page HTML install guide
    #
    # TODO This is currently disabled, because the pandoc-specific
    # buttondown.css doesn't work with the different kind of output
    # makeinfo produces. When we understand better how we want to
    # do generation, decide whether we want multi-page HTML output
    # and how to make it work well.
    #
    # add_custom_command(
    #     OUTPUT ${HTML_OUTPUT_DIR}/install-guide/index.html
    #     COMMAND
    #         ${PANDOC_EXECUTABLE} install-guide.md -o install-guide.texi -s
    #     COMMAND
    #         ${MAKEINFO_EXECUTABLE} install-guide.texi --html -o ${HTML_OUTPUT_DIR}/install-guide --css-ref buttondown.css
    #     DEPENDS
    #         ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
    #         ${HTML_OUTPUT_DIR}/buttondown.css
    #     VERBATIM
    #     )

    # Make the INSTALL file for CPack for the tarball. This gets put
    # into the tarball via the CPack rules in the top-level
    # CMakeLists.txt
    add_custom_command(
        OUTPUT final/INSTALL
        COMMAND
            ${CMAKE_COMMAND} -E make_directory final
        COMMAND
            ${PANDOC_EXECUTABLE} -t plain -o final/INSTALL install-guide.md
        DEPENDS
            ${CMAKE_CURRENT_BINARY_DIR}/install-guide.md
        VERBATIM
        )

    # Add a top-level target for the others to hook onto
    add_custom_target(install-guide
        DEPENDS
            final/INSTALL
            ${HTML_OUTPUT_DIR}/install-guide.html
            ${HTML_OUTPUT_DIR}/install-guide.pdf
        VERBATIM
        )

    # Add a top-level target for the others to hook onto
    add_custom_target(user-guide
        DEPENDS
            ${HTML_OUTPUT_DIR}/user-guide.html
            ${HTML_OUTPUT_DIR}/user-guide.pdf
        VERBATIM
        )
endif()

if(HTML_BUILD_IS_POSSIBLE)
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

    # Add a top-level target for the others to hook onto
    add_custom_target(webpage
        DEPENDS
           ${HTML_OUTPUT_DIR}/index.html
           ${HTML_OUTPUT_DIR}/install-guide.html
           ${HTML_OUTPUT_DIR}/install-guide.pdf
           ${HTML_OUTPUT_DIR}/user-guide.html
           ${HTML_OUTPUT_DIR}/user-guide.pdf
           ${HTML_OUTPUT_DIR}/manual-${PROJECT_VERSION}.pdf
        VERBATIM
        )
    if(BUILD_DOWNLOAD_SECTION)
        add_dependencies(webpage
            ${HTML_OUTPUT_DIR}/gromacs-${PROJECT_VERSION}.tar.gz
            ${HTML_OUTPUT_DIR}/regressiontests-${PROJECT_VERSION}.tar.gz
            )
    endif()
    add_dependencies(webpage webpage-html doc-all)
endif()
