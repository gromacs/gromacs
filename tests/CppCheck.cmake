#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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
#

find_program(CPPCHECK_EXECUTABLE
    NAMES cppcheck
    DOC "cppcheck executable")
mark_as_advanced(CPPCHECK_EXECUTABLE)
set(CPPCHECK_EXTRA_FLAGS "" CACHE STRING "Extra flags to pass to cppcheck")
mark_as_advanced(CPPCHECK_EXTRA_FLAGS)

if (CPPCHECK_EXECUTABLE)
    file(GLOB_RECURSE CPPCHECK_C_FILES RELATIVE ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src/*.c ${CMAKE_SOURCE_DIR}/src/*.cu)
    file(GLOB_RECURSE CPPCHECK_C_FILES_TO_IGNORE RELATIVE ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src/nb_kernel_Elec*.c
        ${CMAKE_SOURCE_DIR}/src/gromacs/linearalgebra/gmx_blas/*.c
        ${CMAKE_SOURCE_DIR}/src/gromacs/linearalgebra/gmx_lapack/*.c
        ${CMAKE_SOURCE_DIR}/src/contrib/*.c
        ${CMAKE_SOURCE_DIR}/src/contrib/*.cu)
    list(REMOVE_ITEM CPPCHECK_C_FILES ${CPPCHECK_C_FILES_TO_IGNORE})
    file(GLOB_RECURSE CPPCHECK_CPP_FILES RELATIVE ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src/*.cpp)

    string(REPLACE ";" "\n" CPPCHECK_C_FILES   "${CPPCHECK_C_FILES}")
    string(REPLACE ";" "\n" CPPCHECK_CPP_FILES "${CPPCHECK_CPP_FILES}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cppcheck-cfiles.txt
        "${CPPCHECK_C_FILES}\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cppcheck-cppfiles.txt
        "${CPPCHECK_CPP_FILES}\n")

    add_custom_target(cppcheck
        COMMAND ${CPPCHECK_EXECUTABLE} ${CPPCHECK_EXTRA_FLAGS} --enable=style
            --file-list=${CMAKE_CURRENT_BINARY_DIR}/cppcheck-cfiles.txt
            -DSIZEOF_LONG_LONG_INT=8 -DSIZEOF_INT=4 -DLINUX
            -I src/gromacs/legacyheaders -I src
            -I src/gromacs/gmxpreprocess
            -I src/programs/mdrun -I src/programs/pdb2gmx
            -I ${CMAKE_BINARY_DIR}/src -I ${CMAKE_BINARY_DIR}/src/gromacs/utility
            --suppress=variableScope
            --suppress=unusedVariable
            --suppress=unreadVariable
            --suppress=unusedStructMember
            --suppress=invalidscanf
            --suppress=sizeofCalculation
            --suppress=missingInclude:src/programs/mdrun/gmx_gpu_utils/gmx_gpu_utils.cu
            --template gcc
            --inline-suppr
            #--xml  2> cppcheck-result.xml
        COMMAND ${CPPCHECK_EXECUTABLE} ${CPPCHECK_EXTRA_FLAGS} --enable=style
            --file-list=${CMAKE_CURRENT_BINARY_DIR}/cppcheck-cppfiles.txt
            -DSIZEOF_LONG_LONG_INT=8 -DSIZEOF_INT=4 -DLINUX
            -I src/gromacs/legacyheaders -I src
            -I ${CMAKE_BINARY_DIR}/src -I ${CMAKE_BINARY_DIR}/src/gromacs/utility
            --suppress=unnecessaryForwardDeclaration
            --suppress=variableScope
            --suppress=missingInclude:src/contrib/openmm_wrapper.cpp
            --suppress=*:src/gromacs/selection/scanner.cpp
            #--xml  2> cppcheck-result-cpp.xml
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running cppcheck"
        VERBATIM)
else()
    add_custom_target(cppcheck
        COMMAND ${CMAKE_COMMAND} -E echo
            "cppcheck was not found by CMake. Rerun CMake specifying CPPCHECK_EXECUTABLE."
        COMMENT "Running cppcheck"
        VERBATIM)
endif()
