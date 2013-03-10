#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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

find_package(Doxygen)
if (DOXYGEN_FOUND)
    # This logic closely follows that found in FindDoxygen.cmake for dot,
    # except that the PATH variable is not cached.
    FIND_PROGRAM(DOXYGEN_MSCGEN_EXECUTABLE
        NAMES mscgen
        DOC "Message sequence chart renderer tool (http://www.mcternan.me.uk/mscgen/)")
    if (DOXYGEN_MSCGEN_EXECUTABLE)
        set(DOXYGEN_MSCGEN_FOUND TRUE)
        get_filename_component(DOXYGEN_MSCGEN_PATH "${DOXYGEN_MSCGEN_EXECUTABLE}" PATH)
    endif (DOXYGEN_MSCGEN_EXECUTABLE)
    mark_as_advanced(DOXYGEN_MSCGEN_EXECUTABLE)
endif (DOXYGEN_FOUND)

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
    if (DOXYGEN_MSCGEN_FOUND)
        set(DOXYGEN_SECTIONS "msc")
    endif (DOXYGEN_MSCGEN_FOUND)
    CONFIGURE_FILE(Doxyfile-common.cmakein Doxyfile-common)
    CONFIGURE_FILE(Doxyfile-full.cmakein Doxyfile-full)
    CONFIGURE_FILE(Doxyfile-lib.cmakein Doxyfile-lib)
    CONFIGURE_FILE(Doxyfile-user.cmakein Doxyfile-user)

    FILE(COPY index.html DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_target(doc-full
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-full
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating full documentation with Doxygen" VERBATIM)
    add_custom_target(doc-lib
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-lib
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating library documentation with Doxygen" VERBATIM)
    add_custom_target(doc-user
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-user
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating user documentation with Doxygen" VERBATIM)
    add_custom_target(doc-all)
    add_dependencies(doc-all doc-full doc-lib doc-user)
endif (DOXYGEN_FOUND)

find_package(PythonInterp)

if (PYTHONINTERP_FOUND AND NOT PYTHON_VERSION_STRING VERSION_LESS "2.6")
    add_custom_target(depgraphs
        ${CMAKE_COMMAND}
        -DSRCDIR=${CMAKE_SOURCE_DIR}
        -DBUILDDIR=${CMAKE_BINARY_DIR}
        -DOUTDIR=${CMAKE_BINARY_DIR}/doxygen/depgraphs
        -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
        -DDOT_EXECUTABLE=${DOXYGEN_DOT_EXECUTABLE}
        -DMODE=GRAPHS
        -P ${CMAKE_SOURCE_DIR}/admin/includedeps.cmake
        COMMENT "Generating include dependency graphs" VERBATIM)
    # TODO: This should be moved away from here and preferably made into
    # a CTest test.
    # TODO: To do the above, issues found by the script should be fixed.
    # A suppression system would also be good to have, since it would be
    # nice to include also legacyheaders in the checks, but those generate
    # a lot of errors about missing documentation.  The script is also
    # currently quite strict, so I don't know whether we really want to
    # enforce all the rules.
    # TODO: Would be nice to not generate the file lists for each of these,
    # but what would be a good shared location for the intermediate files?
    add_custom_target(depcheck
        ${CMAKE_COMMAND}
        -DSRCDIR=${CMAKE_SOURCE_DIR}
        -DBUILDDIR=${CMAKE_BINARY_DIR}
        -DOUTDIR=${CMAKE_BINARY_DIR}/doxygen/depgraphs
        -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
        -DMODE=CHECK
        -P ${CMAKE_SOURCE_DIR}/admin/includedeps.cmake
        COMMENT "Checking include dependencies and documentation" VERBATIM)
endif ()
