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
# images), and content generated from the gmx program for the various
# tools (man and HTML pages). It also provides the "webpage" target,
# that combines all of the above (except man pages in man format) into
# a form suitable for automated deployment to the GROMACS website. It
# also provides the INSTALL file for the tarball.
#
# The webpage is mostly built by Sphinx.  Variable values for Sphinx
# substitutions are configured by CMake (for things like version numbers),
# using gmx_configure_version_file().  This happens during build time instead
# of configure time, because 1) some of the version variables are only
# available during build time, and 2) we don't want to do all the Sphinx setup
# during configuration to save some time when not building the content.
# All the generated values get put into conf-vars.py (generated from
# conf-vars.py.cmakein), which in turn is included by the Sphinx configuration
# file conf.py.

# Even if we aren't going to make the full webpage, set up to put all
# the documentation output in the same place, for convenience
set(HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
file(MAKE_DIRECTORY ${HTML_OUTPUT_DIR})

set(MAN_PAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

find_package(PythonInterp)
find_package(Sphinx 1.2.3 COMPONENTS pygments)

add_subdirectory(doxygen)

option(GMX_BUILD_WEBPAGE "Whether to try to configure to build the GROMACS static webpages" OFF)
mark_as_advanced(GMX_BUILD_WEBPAGE)
option(GMX_BUILD_MANUAL "Whether to try to configure to build the PDF manual" OFF)
mark_as_advanced(GMX_BUILD_MANUAL)
if(GMX_BUILD_MANUAL)
    # Make sure we only do detection of manual-building dependencies
    # when the user opted in for that.
    add_subdirectory(manual)
endif()

if (SPHINX_FOUND)
    # We need to copy all the Sphinx input files to a single directory, and
    # since some of them are (or will be) generated, we just copy everything
    # into the build tree, to this directory.
    set(SPHINX_INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/sphinx-input)
    set(SPHINX_CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/conf.py)
    set(SPHINX_CONFIG_VARS_FILE ${SPHINX_INPUT_DIR}/conf-vars.py)
    set(SPHINX_EXTENSION_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    gmx_configure_version_file(conf-vars.py.cmakein ${SPHINX_CONFIG_VARS_FILE}
        EXTRA_VARS
            SPHINX_EXTENSION_PATH
            GMX_CMAKE_MINIMUM_REQUIRED_VERSION REQUIRED_CUDA_VERSION
            REQUIRED_CUDA_COMPUTE_CAPABILITY REGRESSIONTEST_VERSION
            SOURCE_MD5SUM REGRESSIONTESTS_MD5SUM
        COMMENT "Configuring Sphinx configuration file")
    set(SPHINX_SOURCE_FILES
        index.rst
        download.rst
        install-guide/index.rst
        man/gromacs.7.rst
        user-guide/index.rst
        user-guide/getting-started.rst
        user-guide/flow.rst
        user-guide/cutoff-schemes.rst
        user-guide/mdrun-features.rst
        user-guide/mdrun-performance.rst
        user-guide/mdp-options.rst
        user-guide/file-formats.rst
        user-guide/tools.rst
        user-guide/environment-variables.rst
        user-guide/plotje.gif
        user-guide/xvgr.gif
        conf.py
        links.dat
        )
    set(SPHINX_INPUT_FILES ${SPHINX_CONFIG_VARS_FILE})
    foreach(_file ${SPHINX_SOURCE_FILES})
        add_custom_command(
            OUTPUT ${SPHINX_INPUT_DIR}/${_file}
            COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/${_file}
                ${SPHINX_INPUT_DIR}/${_file}
            DEPENDS
                ${CMAKE_CURRENT_SOURCE_DIR}/${_file}
            COMMENT "Copying Sphinx input file ${_file}"
            VERBATIM)
        list(APPEND SPHINX_INPUT_FILES ${SPHINX_INPUT_DIR}/${_file})
    endforeach()
    gmx_add_custom_output_target(sphinx-input OUTPUT STAMP
        DEPENDS ${SPHINX_INPUT_FILES})
    gmx_add_custom_output_target(sphinx-programs OUTPUT STAMP
        COMMAND ${CMAKE_COMMAND} -E make_directory programs
        COMMAND gmx -quiet help -export rst
        DEPENDS gmx
        WORKING_DIRECTORY ${SPHINX_INPUT_DIR}
        COMMENT "Generating reStructuredText help")
    # This dependency ensures that the directories exist before the
    # executable tries to write things there.
    add_dependencies(sphinx-programs sphinx-input)

    # Make the INSTALL file for CPack for the tarball. This gets put
    # into the tarball via the CPack rules below, which requires that
    # the INSTALL file is in a separate directory by itself.
    set(TEXT_INSTALL_GUIDE_OUTPUT_DIR "install-guide/text")
    add_custom_target(install-guide
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b text
            -w sphinx-install.log
            -d ${CMAKE_CURRENT_BINARY_DIR}/install-guide/_doctrees
            -c ${SPHINX_INPUT_DIR}
            "${SPHINX_INPUT_DIR}/install-guide"
            "${TEXT_INSTALL_GUIDE_OUTPUT_DIR}"
        COMMAND
            ${CMAKE_COMMAND} -E rename
            ${TEXT_INSTALL_GUIDE_OUTPUT_DIR}/index.txt
            ${TEXT_INSTALL_GUIDE_OUTPUT_DIR}/INSTALL
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building INSTALL with Sphinx"
        VERBATIM
        )
    add_dependencies(install-guide sphinx-input)
    gmx_cpack_add_generated_source_directory(install-guide/text DESTINATION /)

    # For a real build of the webpage, the md5sum of the tarballs must
    # already be known, and so we may as well require that the real
    # build of the webpage take place from cmake run from the unpacked
    # tarball. Then, the *_MD5SUM and *_TARBALL variables will be able
    # to be set on the cmake command line (e.g. by a Jenkins job
    # configuration), and we can require that they are set. For local
    # building of the webpages (e.g. for debugging), those variables
    # can be left unset, and if so, the download section will not be
    # constructed.
    if(SOURCE_IS_SOURCE_DISTRIBUTION AND GMX_BUILD_WEBPAGE)
        foreach(VAR SOURCE_MD5SUM REGRESSIONTESTS_MD5SUM SOURCE_TARBALL REGRESSIONTESTS_TARBALL)
            if(NOT DEFINED ${VAR})
                message(FATAL_ERROR "The build of the webpage requires that ${VAR} is set in the cmake cache, e.g. on the CMake command line")
            endif()
        endforeach()
        set(BUILD_DOWNLOAD_SECTION on)
    else()
        if (SOURCE_TARBALL AND SOURCE_MD5SUM AND
                REGRESSIONTESTS_TARBALL AND REGRESSIONTESTS_MD5SUM)
            set(BUILD_DOWNLOAD_SECTION on)
        else()
            set(BUILD_DOWNLOAD_SECTION off)
        endif()
    endif()

    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

    add_custom_target(webpage-sphinx
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${SPHINX_INPUT_DIR}/_static
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b html
            -w sphinx-html.log
            -d "${SPHINX_CACHE_DIR}"
            "${SPHINX_INPUT_DIR}"
            "${HTML_OUTPUT_DIR}"
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building HTML documentation with Sphinx"
        VERBATIM
        )
    add_dependencies(webpage-sphinx sphinx-input sphinx-programs)

    add_custom_target(man
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b man
            -w sphinx-man.log
            -d ${SPHINX_CACHE_DIR}
            -t do_man
            ${SPHINX_INPUT_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/man
        COMMENT "Building man pages with Sphinx"
        VERBATIM)
    add_dependencies(man sphinx-input sphinx-programs)
    if (NOT SOURCE_IS_SOURCE_DISTRIBUTION OR GMX_BUILD_HELP)
        set(MAN_PAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/man)
    endif()
endif()
if (SOURCE_IS_SOURCE_DISTRIBUTION OR (GMX_BUILD_HELP AND SPHINX_FOUND))
    install(DIRECTORY ${MAN_PAGE_DIR}/
        DESTINATION ${MAN_INSTALL_DIR}/man1
        COMPONENT man OPTIONAL
        FILES_MATCHING PATTERN "*.1")
    install(DIRECTORY ${MAN_PAGE_DIR}/
        DESTINATION ${MAN_INSTALL_DIR}/man7
        COMPONENT man OPTIONAL
        FILES_MATCHING PATTERN "*.7")
endif()
gmx_cpack_add_generated_source_directory(man)

set(HTML_BUILD_IS_POSSIBLE OFF)
# We can only configure to build the webpage if the user asked for it,
# the build is outside of the source dir, and all the components can
# be built. There's no intrinsic need to be talkative if we fail -
# most people never need to know, unless they've asked for the webpage
# build.
if(GMX_BUILD_WEBPAGE)
    set(HTML_BUILD_IS_POSSIBLE ON)
    # Next, turn it off in any of the preconditions are unsatisified

    if(NOT MANUAL_BUILD_IS_POSSIBLE)
        set(HTML_BUILD_IS_POSSIBLE OFF)
        message(STATUS "Cannot build webpage without being able to build the reference PDF manual")
    endif()

    if(NOT PYTHON_EXECUTABLE)
        set(HTML_BUILD_IS_POSSIBLE OFF)
        message(STATUS "Cannot build webpage without python")
    endif()

    if(NOT SPHINX_FOUND)
        set(HTML_BUILD_IS_POSSIBLE OFF)
        message(STATUS "Cannot build webpage without sphinx")
    endif()

    if(NOT DOXYGEN_EXECUTABLE)
        set(HTML_BUILD_IS_POSSIBLE OFF)
        message(STATUS "Cannot build webpage without doxygen")
    endif()

    if(NOT DOXYGEN_MSCGEN_EXECUTABLE)
        set(HTML_BUILD_IS_POSSIBLE OFF)
        message(STATUS "Cannot build webpage without mscgen")
    endif()
endif()

# TODO Get rid of the remaining content in old-html
add_custom_target(webpage-old-images
    ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/old-html ${HTML_OUTPUT_DIR}
    COMMENT "Copying flowchart images"
    VERBATIM
    )

if(HTML_BUILD_IS_POSSIBLE)
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

    # The Doxygen configuration in doxygen/Doxyfile-common.cmakein
    # makes all the Doxygen output directly in
    # ${HTML_OUTPUT_DIR}/doxygen (and makes the directory if it needs
    # to).

    # Add a top-level target for the others to hook onto
    add_custom_target(webpage
        DEPENDS
           ${HTML_OUTPUT_DIR}/manual-${GMX_VERSION_STRING}.pdf
           ${extra_webpage_dependencies}
        VERBATIM
        )
    add_dependencies(webpage webpage-sphinx webpage-old-images doxygen-all)
else()
    add_custom_target(webpage
        COMMAND ${CMAKE_COMMAND} -E echo "Cannot build webpage because of missing requirements. Check cmake status output for reasons"
        COMMENT "Webpage build disabled"
        )
endif()
