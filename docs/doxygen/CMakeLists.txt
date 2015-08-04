#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

include(gmxCustomCommandUtilities)
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

# Note that if you relocate this file, you also need to change gmxtree.py.
# TODO: Try to find a solution where this would be hard-coded here, completely
# separate from the logic that generates it.
set(INSTALLED_HEADERS_FILE ${PROJECT_BINARY_DIR}/src/gromacs/installed-headers.txt)

########################################################################
# Doxygen configuration
#
if (DOXYGEN_FOUND)
    # Generate Doxygen input configuration files
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

    # Add a stamp that gets updated whenever a source file that contributes to
    # Doxygen is updated.  All Doxygen targets can then depend on this stamp,
    # instead of each specifying this massive set of dependencies.
    # TODO: The set of files that affects Doxygen is not exactly the same as
    # that which affects the dependency graphs, and neither is captured
    # perfectly with this set, but this can be improved if it becomes a real
    # issue (the differences are mainly in the generated group kernel files,
    # but there are some minor other differences as well).
    if (SOURCE_IS_GIT_REPOSITORY AND GIT_FOUND)
        gmx_get_files_with_gitattribute(doxygen_sources gmx-doxygen)
        # Add also common configuration files that influence the generation.
        set(doxygen_deps
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-common
            ${CMAKE_CURRENT_BINARY_DIR}/RunDoxygen.cmake)
        if (GMX_COMPACT_DOXYGEN)
            list(APPEND doxygen_deps ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile-compact)
        endif()
        gmx_add_custom_output_target(doxygen-source-timestamp
            OUTPUT STAMP
            DEPENDS_FILE_LIST doxygen_sources
            DEPENDS ${doxygen_deps}
            COMMENT "Updating Doxygen source timestamp")
    else()
        # If there is no git, we just update the stamp every time, and the
        # builds are triggered every time.  This should be rare enough, but
        # this case still needs to be supported for the Jenkins job that builds
        # the documentation for a release from the release tarball.
        gmx_add_custom_output_target(doxygen-source-timestamp
            RUN_ALWAYS OUTPUT STAMP
            COMMENT "Updating Doxygen source timestamp")
    endif()

    # Helper function to create a Doxygen build target
    function (add_doxygen_target TARGET TYPE COMMENT)
        set(_deps doxygen-source-timestamp doxygen-version
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-${TYPE})
        set(_fast_arg "")
        if ("${ARGN}" STREQUAL "DEPGRAPHS")
            set(_fast_arg "ADD_FAST_TARGET")
            if (TARGET dep-graphs-dot)
                list(APPEND _deps dep-graphs-dot)
            endif()
        endif()
        gmx_add_custom_output_target(${TARGET} OUTPUT STAMP ${_fast_arg}
            # Ensure the directory exists to avoid spurious warnings
            COMMAND ${CMAKE_COMMAND} -E make_directory ${DEPGRAPH_DIR}
            COMMAND ${CMAKE_COMMAND} -DDOCTYPE=${TYPE} -P RunDoxygen.cmake
            DEPENDS ${_deps}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "${COMMENT}")
        if (_fast_arg)
            add_dependencies(${TARGET}-fast doxygen-version)
        endif()
    endfunction()

    # The XML target is necessary for the dependency graphs
    add_doxygen_target(doxygen-xml xml "Extracting Doxygen documentation to XML")

    # Add targets for generating the dependency graphs; they are necessary for
    # the full and lib variants of the documentation.
    if (USE_PYTHON_SCRIPTS)
        # Add target that generates .dot files for the dependency graphs.
        # The output also depends on the Doxygen source files, but the
        # dependency on doxygen-xml takes care of that transitively.
        gmx_add_custom_output_target(dep-graphs-dot ADD_FAST_TARGET
            OUTPUT STAMP
            COMMAND ${PYTHON_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/graphbuilder.py
                -S ${CMAKE_SOURCE_DIR} -B ${CMAKE_BINARY_DIR}
                --ignore-cycles ${CMAKE_CURRENT_SOURCE_DIR}/cycle-suppressions.txt
                -o ${DEPGRAPH_DIR}
            DEPENDS doxygen-xml ${INSTALLED_HEADER_FILE}
                ${CMAKE_CURRENT_SOURCE_DIR}/doxygenxml.py
                ${CMAKE_CURRENT_SOURCE_DIR}/gmxtree.py
                ${CMAKE_CURRENT_SOURCE_DIR}/graphbuilder.py
                ${CMAKE_CURRENT_SOURCE_DIR}/cycle-suppressions.txt
            COMMENT "Generating include dependency graphs for dot")

        # Add separate targets for converting the .dot files to .png.
        # These are not needed by the Doxygen targets, but are useful for
        # investigating just the graphs.
        set(dep_graphs_command_dot
            ${CMAKE_COMMAND} -DGRAPHDIR=${DEPGRAPH_DIR}
            -DDOT_EXECUTABLE=${DOXYGEN_DOT_EXECUTABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/generateGraphs.cmake)
        add_custom_target(dep-graphs      COMMAND ${dep_graphs_command_dot}
            COMMENT "Generating PNG include dependency graphs" VERBATIM)
        add_custom_target(dep-graphs-fast COMMAND ${dep_graphs_command_dot}
            COMMENT "Generating PNG include dependency graphs" VERBATIM)
        add_dependencies(dep-graphs      dep-graphs-dot)
        add_dependencies(dep-graphs-fast dep-graphs-dot-fast)
    endif()

    # Create the actual targets to build the documentation.
    add_doxygen_target(doxygen-full full "Generating full documentation with Doxygen" DEPGRAPHS)
    add_doxygen_target(doxygen-lib  lib  "Generating library documentation with Doxygen" DEPGRAPHS)
    add_doxygen_target(doxygen-user user "Generating public API documentation with Doxygen")
    # Convenience targets to build all the documentation flavors.
    add_custom_target(doxygen-all)
    add_custom_target(doxygen-all-fast)
    add_dependencies(doxygen-all      doxygen-full doxygen-lib doxygen-user)
    add_dependencies(doxygen-all-fast doxygen-full-fast doxygen-lib-fast doxygen-user)

    # Finally, create the check-source target.
    if (USE_PYTHON_SCRIPTS)
        # TODO: Consider whether this is the best location for this code,
        # since not all of it is Doxygen-specific (but nearly all of it
        # relies on the Doxygen XML documentation).
        # The output .log file currently needs to be here, since Jenkins
        # expects that.
        set(check_source_command
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check-source.py
            -S ${CMAKE_SOURCE_DIR} -B ${CMAKE_BINARY_DIR}
            -l ${CMAKE_CURRENT_BINARY_DIR}/check-source.log
            --ignore ${CMAKE_CURRENT_SOURCE_DIR}/suppressions.txt
            --ignore-cycles ${CMAKE_CURRENT_SOURCE_DIR}/cycle-suppressions.txt)
        add_custom_target(check-source      COMMAND ${check_source_command}
            COMMENT "Checking source code for various issues" VERBATIM)
        add_custom_target(check-source-fast COMMAND ${check_source_command}
            COMMENT "Checking source code for various issues" VERBATIM)
        add_dependencies(check-source       doxygen-xml)
    endif()
else()
    function (add_disabled_doxygen_target TARGET)
        add_custom_target(${TARGET}
            COMMAND ${CMAKE_COMMAND} -E echo
                "Cannot build/run target ${TARGET} because Doxygen is not available"
            COMMENT "Doxygen not available" VERBATIM)
    endfunction()
    add_disabled_doxygen_target(doxygen-all)
    add_disabled_doxygen_target(doxygen-all-fast)
    add_disabled_doxygen_target(doxygen-xml)
    add_disabled_doxygen_target(doxygen-user)
    add_disabled_doxygen_target(doxygen-lib)
    add_disabled_doxygen_target(doxygen-lib-fast)
    add_disabled_doxygen_target(doxygen-full)
    add_disabled_doxygen_target(doxygen-full-fast)
    add_disabled_doxygen_target(check-source)
    add_disabled_doxygen_target(check-source-fast)
    add_disabled_doxygen_target(dep-graphs)
    add_disabled_doxygen_target(dep-graphs-fast)
endif()
