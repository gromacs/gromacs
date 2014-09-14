#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

include(gmxOptionUtilities)

find_package(Doxygen QUIET)
if (DOXYGEN_FOUND)
    # This logic closely follows that found in FindDoxygen.cmake for dot,
    # except that the PATH variable is not cached.
    find_program(DOXYGEN_MSCGEN_EXECUTABLE
        NAMES mscgen
        DOC "Message sequence chart renderer tool (http://www.mcternan.me.uk/mscgen/)")
    if (DOXYGEN_MSCGEN_EXECUTABLE)
        set(DOXYGEN_MSCGEN_FOUND TRUE)
        get_filename_component(DOXYGEN_MSCGEN_PATH "${DOXYGEN_MSCGEN_EXECUTABLE}" PATH)
    endif()
    mark_as_advanced(DOXYGEN_MSCGEN_EXECUTABLE)
endif()

gmx_dependent_option(
    GMX_COMPACT_DOXYGEN
    "Skip parts of Doxygen documentation (in particular, many graphs) to save disk space and time"
    OFF
    DOXYGEN_FOUND)
mark_as_advanced(GMX_COMPACT_DOXYGEN)

find_package(PythonInterp)
set(USE_PYTHON_SCRIPTS OFF)
if (PYTHONINTERP_FOUND AND NOT PYTHON_VERSION_STRING VERSION_LESS "2.6")
    set(USE_PYTHON_SCRIPTS ON)
endif()

# Create a stamp file whenever cmake is run, as the list of installed files may
# have changed.
set(CONFIGURE_STAMP_FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake-stamp.txt)
file(WRITE ${CONFIGURE_STAMP_FILE} "Stamp for configuration")
# Note that if you relocate this file, you also need to change gmxtree.py.
set(INSTALLED_HEADERS_FILE ${CMAKE_CURRENT_BINARY_DIR}/installed-headers.txt)
add_custom_command(OUTPUT ${INSTALLED_HEADERS_FILE}
    COMMAND ${CMAKE_COMMAND}
        -D SRCDIR=${CMAKE_SOURCE_DIR}
        -D BUILDDIR=${CMAKE_BINARY_DIR}
        -D OUTFILE=${INSTALLED_HEADERS_FILE}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/getInstalledHeaders.cmake
    DEPENDS ${CONFIGURE_STAMP_FILE}
    COMMENT "Generating installed headers list" VERBATIM)
add_custom_target(find-installed-headers DEPENDS ${INSTALLED_HEADERS_FILE})

########################################################################
# Doxygen configuration
########################################################################
if (DOXYGEN_FOUND)
    file(GLOB NB_KERNEL_DIRS
         ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/nonbonded/nb_kernel_*/)
    list(REMOVE_ITEM NB_KERNEL_DIRS
         ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/nonbonded/nb_kernel_c)
    foreach (dir ${NB_KERNEL_DIRS})
        SET(NB_KERNEL_DIRS_TO_IGNORE_IN_DOXYGEN
            "${NB_KERNEL_DIRS_TO_IGNORE_IN_DOXYGEN} \\\n                         ${dir}")
    endforeach()
    set(DEPGRAPH_DIR ${CMAKE_CURRENT_BINARY_DIR}/depgraphs)
    set(DOXYGEN_SECTIONS "")
    set(DOXYGEN_EXTRA_SETTINGS "")
    if (GMX_COMPACT_DOXYGEN)
        set(DOXYGEN_EXTRA_SETTINGS "@INCLUDE   = ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile-compact")
    endif()
    configure_file(Doxyfile-common.cmakein Doxyfile-common)
    configure_file(Doxyfile-full.cmakein Doxyfile-full)
    configure_file(Doxyfile-lib.cmakein Doxyfile-lib)
    configure_file(Doxyfile-user.cmakein Doxyfile-user)
    configure_file(Doxyfile-xml.cmakein Doxyfile-xml)

    configure_file(RunDoxygen.cmake.cmakein RunDoxygen.cmake @ONLY)
    gmx_configure_version_file(Doxyfile-version.cmakein Doxyfile-version
                               TARGET doxygen-version)

    set(SOURCE_STAMP_FILE ${CMAKE_CURRENT_BINARY_DIR}/source-stamp.txt)
    set(source_stamp_depends "")
    if (SOURCE_IS_GIT_REPOSITORY AND GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} ls-files
            COMMAND ${GIT_EXECUTABLE} check-attr gmx-doxygen --stdin
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE doxygen_dep_attr)
        string(REGEX MATCHALL "[^\n]*: gmx-doxygen: set\n"
               doxygen_deps "${doxygen_dep_attr}")
        string(REGEX REPLACE "([^;]*): gmx-doxygen: set\n" "${PROJECT_SOURCE_DIR}/\\1"
               doxygen_deps "${doxygen_deps}")
        list(APPEND doxygen_deps ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml)
        list(APPEND doxygen_deps ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-common)
        list(APPEND doxygen_deps ${CMAKE_CURRENT_BINARY_DIR}/RunDoxygen.cmake)
        if (GMX_COMPACT_DOXYGEN)
            list(APPEND doxygen_deps ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-compact)
        endif()
        add_custom_command(OUTPUT ${SOURCE_STAMP_FILE}
            COMMAND ${CMAKE_COMMAND} -E touch ${SOURCE_STAMP_FILE}
            DEPENDS ${doxygen_deps}
            COMMENT "Setting Doxygen source stamp" VERBATIM)
        add_custom_target(doxygen-source-stamp DEPENDS ${SOURCE_STAMP_FILE})
        if (CMAKE_GENERATOR STREQUAL "Ninja")
            set(source_stamp_depends ${SOURCE_STAMP_FILE})
        endif()
    else()
        add_custom_target(doxygen-source-stamp
            COMMAND ${CMAKE_COMMAND} -E touch ${SOURCE_STAMP_FILE}
            VERBATIM)
    endif()

    function (add_doxygen_target TARGET TYPE COMMENT)
        set(_stamp_file doxygen-${TYPE}-stamp.txt)
        set(_deps ${source_stamp_depends} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-${TYPE})
        if ("${ARGN}" STREQUAL "DEPGRAPHS")
            list(APPEND _deps dep-graphs-dot)
            if (CMAKE_GENERATOR STREQUAL "Ninja")
                list(APPEND _deps ${dep_graphs_stamp})
            endif()
        endif()
        add_custom_command(OUTPUT ${_stamp_file}
            # Ensure the directory exists to avoid spurious warnings
            COMMAND ${CMAKE_COMMAND} -E make_directory ${DEPGRAPH_DIR}
            COMMAND ${CMAKE_COMMAND} -DDOCTYPE=${TYPE} -P RunDoxygen.cmake
            COMMAND ${CMAKE_COMMAND} -E touch ${_stamp_file}
            DEPENDS ${_deps}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "${COMMENT}" VERBATIM)
        add_custom_target(${TARGET} DEPENDS ${_stamp_file})
        add_dependencies(${TARGET} doxygen-version doxygen-source-stamp)
    endfunction()

    add_doxygen_target(doc-xml       xml  "Extracting Doxygen documentation to XML")

    set(dep_graphs_stamp "")
    if (USE_PYTHON_SCRIPTS)
        # TODO: Consider whether this is the best name and location for this
        # code, since not all of it is Doxygen-specific (but nearly all of it
        # relies on the Doxygen XML documentation).
        set(doc_check_command
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-check.py
            -S ${CMAKE_SOURCE_DIR} -B ${CMAKE_BINARY_DIR}
            -l ${CMAKE_CURRENT_BINARY_DIR}/doxygen-check.log
            --exitcode
            --ignore ${CMAKE_CURRENT_SOURCE_DIR}/suppressions.txt
            --ignore-cycles ${CMAKE_CURRENT_SOURCE_DIR}/cycle-suppressions.txt)
        add_custom_target(doc-check COMMAND ${doc_check_command}
            COMMENT "Checking Doxygen documentation" VERBATIM)
        add_dependencies(doc-check doc-xml find-installed-headers)
        add_custom_target(doc-check-fast COMMAND ${doc_check_command}
            COMMENT "Checking Doxygen documentation" VERBATIM)
        add_dependencies(doc-check-fast find-installed-headers)

        set(dep_graphs_stamp ${DEPGRAPH_DIR}/dot-stamp.txt)
        set(dep_graphs_deps
            ${CMAKE_CURRENT_SOURCE_DIR}/doxygenxml.py
            ${CMAKE_CURRENT_SOURCE_DIR}/gmxtree.py
            ${CMAKE_CURRENT_SOURCE_DIR}/graphbuilder.py
            ${CMAKE_CURRENT_SOURCE_DIR}/cycle-suppressions.txt
            doc-xml find-installed-headers)
        if (CMAKE_GENERATOR STREQUAL "Ninja")
            list(APPEND dep_graphs_deps
                 ${CMAKE_CURRENT_BINARY_DIR}/doxygen-xml-stamp.txt
                 ${INSTALLED_HEADERS_FILE})
        endif()
        set(dep_graphs_command_python
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/graphbuilder.py
            -S ${CMAKE_SOURCE_DIR} -B ${CMAKE_BINARY_DIR}
            --ignore-cycles ${CMAKE_CURRENT_SOURCE_DIR}/cycle-suppressions.txt
            -o ${DEPGRAPH_DIR})
        set(dep_graphs_command_dot
            ${CMAKE_COMMAND} -DGRAPHDIR=${DEPGRAPH_DIR}
            -DDOT_EXECUTABLE=${DOXYGEN_DOT_EXECUTABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/generateGraphs.cmake)
        add_custom_command(OUTPUT ${dep_graphs_stamp}
            COMMAND ${dep_graphs_command_python}
            COMMAND ${CMAKE_COMMAND} -E touch ${dep_graphs_stamp}
            DEPENDS ${dep_graphs_deps}
            COMMENT "Generating include dependency graphs for dot" VERBATIM)
        add_custom_target(dep-graphs-dot DEPENDS ${dep_graphs_stamp})
        add_custom_target(dep-graphs
            COMMAND ${dep_graphs_command_dot}
            COMMENT "Generating PNG include dependency graphs" VERBATIM)
        add_dependencies(dep-graphs dep-graphs-dot)
        add_custom_target(dep-graphs-dot-fast
            COMMAND ${dep_graphs_command_python}
            COMMENT "Generating include dependency graphs for dot" VERBATIM)
        add_custom_target(dep-graphs-fast
            COMMAND ${dep_graphs_command_dot}
            COMMENT "Generating PNG include dependency graphs" VERBATIM)
        add_dependencies(dep-graphs-dot-fast find-installed-headers)
        add_dependencies(dep-graphs-fast dep-graphs-dot-fast)
    endif()

    # The targets ending with -fast do the same thing as the target without the
    # suffix, but assume that time-consuming dependencies have already been
    # built, making it faster and more convenient to test a single part of the
    # system.
    add_doxygen_target(doc-full      full "Generating full documentation with Doxygen" DEPGRAPHS)
    add_doxygen_target(doc-full-fast full "Generating full documentation with Doxygen")
    add_doxygen_target(doc-lib       lib  "Generating library documentation with Doxygen" DEPGRAPHS)
    add_doxygen_target(doc-lib-fast  lib  "Generating library documentation with Doxygen")
    add_doxygen_target(doc-user      user "Generating public API documentation with Doxygen")
    add_custom_target(doc-all)
    add_custom_target(doc-all-fast)
    add_dependencies(doc-all doc-full doc-lib doc-user)
    add_dependencies(doc-all-fast doc-full-fast doc-lib-fast doc-user)
endif()
