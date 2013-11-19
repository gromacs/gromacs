#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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

function (generate_module_file_list SRCDIR OUTFILE MODE)
    set(_module_list
        analysisdata commandline fft fileio linearalgebra onlinehelp options
        selection timing trajectoryanalysis utility)
    if (MODE STREQUAL "CHECK")
        list(APPEND _module_list gmxana gmxlib gmxpreprocess legacyheaders mdlib)
    endif()
    set(PATH_LIST)
    foreach (MODULE ${_module_list})
        list(APPEND PATH_LIST "${SRCDIR}/src/gromacs/${MODULE}/*.cpp")
        if (MODE STREQUAL "GRAPHS")
            list(APPEND PATH_LIST "${SRCDIR}/src/gromacs/${MODULE}/*.c")
        endif()
        list(APPEND PATH_LIST "${SRCDIR}/src/gromacs/${MODULE}/*.h")
    endforeach ()
    list(APPEND PATH_LIST "${SRCDIR}/src/testutils/*.cpp")
    list(APPEND PATH_LIST "${SRCDIR}/src/testutils/*.h")
    set(FILE_LIST)
    foreach (PATH_EXPR ${PATH_LIST})
        file(GLOB_RECURSE FOUND_FILES ${PATH_EXPR})
        list(APPEND FILE_LIST ${FOUND_FILES})
    endforeach ()
    string(REPLACE ";" "\n" FILE_LIST "${FILE_LIST}")
    file(WRITE ${OUTFILE} "${FILE_LIST}")
endfunction ()

function (generate_installed_file_list SRCDIR BUILDDIR OUTFILE)
    file(GLOB_RECURSE INSTALL_FILE_LIST "${BUILDDIR}/cmake_install.cmake")
    set(MATCH_REGEX "${SRCDIR}/.*\\.h")
    set(HEADER_LIST)
    foreach (INSTALL_FILE ${INSTALL_FILE_LIST})
        file(STRINGS ${INSTALL_FILE} HEADER_LINES REGEX "${MATCH_REGEX}")
        foreach (HEADER_LINE ${HEADER_LINES})
            string (REGEX MATCH "${MATCH_REGEX}" HEADER "${HEADER_LINE}")
            list(APPEND HEADER_LIST "${HEADER}")
        endforeach ()
    endforeach ()
    string(REPLACE ";" "\n" HEADER_LIST "${HEADER_LIST}")
    file(WRITE ${OUTFILE} "${HEADER_LIST}")
endfunction ()

if (NOT DEFINED SRCDIR OR NOT DEFINED BUILDDIR OR NOT DEFINED OUTDIR)
    message(FATAL_ERROR "Required input variable (SRCDIR, BUILDDIR, OUTDIR) not set")
endif ()

if (NOT DEFINED PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE python)
endif ()

if (NOT DEFINED MODE)
    set(MODE "CHECK")
endif ()

if (MODE STREQUAL "CHECK")
    set(GRAPHOPTIONS --check)
elseif (MODE STREQUAL "CHECKDOC")
    set(GRAPHOPTIONS --check --check-doc --warn-undoc)
elseif (MODE STREQUAL "GRAPHS")
    set(GRAPHOPTIONS
        --module-graph module-deps.dot --module-file-graphs
        -o ${OUTDIR})
else ()
    message(FATAL_ERROR "Unknown mode ${MODE}")
endif ()

file(MAKE_DIRECTORY ${OUTDIR})
generate_module_file_list(${SRCDIR} ${OUTDIR}/module-files.txt ${MODE})
generate_installed_file_list(${SRCDIR} ${BUILDDIR} ${OUTDIR}/installed-headers.txt)
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SRCDIR}/admin/includedeps.py
                        -f ${OUTDIR}/module-files.txt
                        --installed ${OUTDIR}/installed-headers.txt
                        -R ${SRCDIR}/src -R ${BUILDDIR}/src
                        -I ${SRCDIR}/src/gromacs/legacyheaders
                        -I ${BUILDDIR}/src/gromacs/utility
                        ${GRAPHOPTIONS})

if (MODE STREQUAL "GRAPHS" AND DOT_EXECUTABLE)
    file(GLOB DOT_INPUT_FILES ${OUTDIR}/*.dot)
    execute_process(COMMAND ${DOT_EXECUTABLE} -Tpng -O ${DOT_INPUT_FILES})
endif ()
