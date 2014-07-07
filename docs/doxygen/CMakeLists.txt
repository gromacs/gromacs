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
    FIND_PROGRAM(DOXYGEN_MSCGEN_EXECUTABLE
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

add_custom_target(find-installed-headers
    COMMAND ${CMAKE_COMMAND}
        -D SRCDIR=${CMAKE_SOURCE_DIR}
        -D BUILDDIR=${CMAKE_BINARY_DIR}
        -D OUTFILE=${CMAKE_CURRENT_BINARY_DIR}/installed-headers.txt
        -P ${CMAKE_CURRENT_SOURCE_DIR}/getInstalledHeaders.cmake
    COMMENT "Generating installed headers list" VERBATIM)

########################################################################
# Doxygen configuration
########################################################################
if (DOXYGEN_FOUND)
    FILE(GLOB NB_KERNEL_DIRS
         ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/nonbonded/nb_kernel_*/)
    LIST(REMOVE_ITEM NB_KERNEL_DIRS
         ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/nonbonded/nb_kernel_c)
    FOREACH(NB_KERNEL_DIR ${NB_KERNEL_DIRS})
        SET(NB_KERNEL_DIRS_TO_IGNORE_IN_DOXYGEN
            "${NB_KERNEL_DIRS_TO_IGNORE_IN_DOXYGEN} \\\n                         ${NB_KERNEL_DIR}")
    ENDFOREACH(NB_KERNEL_DIR)
    set(DOXYGEN_SECTIONS "")
    CONFIGURE_FILE(Doxyfile-common.cmakein Doxyfile-common)
    CONFIGURE_FILE(Doxyfile-full.cmakein Doxyfile-full)
    CONFIGURE_FILE(Doxyfile-lib.cmakein Doxyfile-lib)
    CONFIGURE_FILE(Doxyfile-user.cmakein Doxyfile-user)
    CONFIGURE_FILE(Doxyfile-xml.cmakein Doxyfile-xml)

    if (GMX_COMPACT_DOXYGEN)
        FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-common
             "@INCLUDE   = ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile-compact\n")
    endif()

    configure_file(RunDoxygen.cmake.cmakein RunDoxygen.cmake @ONLY)
    add_custom_target(doc-full
        ${CMAKE_COMMAND} -DDOCTYPE=full -P RunDoxygen.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating full documentation with Doxygen" VERBATIM)
    add_custom_target(doc-lib
        ${CMAKE_COMMAND} -DDOCTYPE=lib -P RunDoxygen.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating library documentation with Doxygen" VERBATIM)
    add_custom_target(doc-user
        ${CMAKE_COMMAND} -DDOCTYPE=user -P RunDoxygen.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating public API documentation with Doxygen" VERBATIM)
    add_custom_target(doc-xml
        ${CMAKE_COMMAND} -DDOCTYPE=xml -P RunDoxygen.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Extracting Doxygen documentation to XML" VERBATIM)
    add_custom_target(doc-all)
    add_dependencies(doc-all doc-full doc-lib doc-user)

    if (GMX_GIT_VERSION_INFO)
        add_custom_target(doxygen-version
                COMMAND ${CMAKE_COMMAND}
                    -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                    -D PROJECT_VERSION="${PROJECT_VERSION}"
                    -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                    -D VERSION_CMAKEIN="${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile-version.cmakein"
                    -D VERSION_OUT="${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-version"
                    -D VERSION_NO_REMOTE_HASH=
                    -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile-version.cmakein
                COMMENT "Generating version information for Doxygen")
        add_dependencies(doc-full doxygen-version)
        add_dependencies(doc-lib doxygen-version)
        add_dependencies(doc-user doxygen-version)
        add_dependencies(doc-xml doxygen-version)
    else()
        set(GMX_PROJECT_VERSION_STR ${PROJECT_VERSION})
        configure_file(Doxyfile-version.cmakein Doxyfile-version)
    endif()

    if (USE_PYTHON_SCRIPTS)
        # TODO: Consider whether this is the best name and location for this
        # code, since not all of it is Doxygen-specific (but nearly all of it
        # relies on the Doxygen XML documentation).
        set(doc_check_command
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-check.py
            -S ${CMAKE_SOURCE_DIR} -B ${CMAKE_BINARY_DIR}
            --installed ${CMAKE_CURRENT_BINARY_DIR}/installed-headers.txt
            -l ${CMAKE_CURRENT_BINARY_DIR}/doxygen-check.log
            --ignore ${CMAKE_CURRENT_SOURCE_DIR}/suppressions.txt)
        add_custom_target(doc-check COMMAND ${doc_check_command}
            COMMENT "Checking Doxygen documentation" VERBATIM)
        add_dependencies(doc-check doc-xml find-installed-headers)

        set(graphdir ${CMAKE_CURRENT_BINARY_DIR}/depgraphs)
        set(dep_graphs_command_python
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/graphbuilder.py
            -S ${CMAKE_SOURCE_DIR} -B ${CMAKE_BINARY_DIR}
            --installed ${CMAKE_CURRENT_BINARY_DIR}/installed-headers.txt
            -o ${CMAKE_CURRENT_BINARY_DIR}/depgraphs)
        set(dep_graphs_command_dot
            ${CMAKE_COMMAND} -DGRAPHDIR=${graphdir}
            -DDOT_EXECUTABLE=${DOXYGEN_DOT_EXECUTABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/generateGraphs.cmake)
        add_custom_target(dep-graphs
            COMMAND ${dep_graphs_command_python}
            COMMAND ${dep_graphs_command_dot}
            COMMENT "Generating include dependency graphs" VERBATIM)
        add_dependencies(dep-graphs doc-xml find-installed-headers)

        # These targets are the same as above, but they don't rerun the
        # dependencies each time, making it faster and more convenient for
        # testing.
        add_custom_target(doc-check-fast COMMAND ${doc_check_command}
            COMMENT "Checking Doxygen documentation" VERBATIM)
        add_custom_target(dep-graphs-fast
            COMMAND ${dep_graphs_command_python}
            COMMAND ${dep_graphs_command_dot}
            COMMENT "Generating include dependency graphs" VERBATIM)
    endif()
endif()
