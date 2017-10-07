#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

find_program(CPPCHECK_EXECUTABLE
    NAMES cppcheck
    DOC "cppcheck executable")
mark_as_advanced(CPPCHECK_EXECUTABLE)
option(CPPCHECK_XML_OUTPUT "Whether to produce XML output from cppcheck (mainly for Jenkins)"
       OFF)
mark_as_advanced(CPPCHECK_XML_OUTPUT)

# Depends on stderr redirection and cat
if (CPPCHECK_EXECUTABLE AND UNIX)
    # Produce the list of files to check
    file(GLOB_RECURSE _inputfiles RELATIVE ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src/*.c
        ${CMAKE_SOURCE_DIR}/src/*.cpp
        ${CMAKE_SOURCE_DIR}/src/*.cu)
    file(GLOB_RECURSE _files_to_ignore RELATIVE ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src/nb_kernel_Elec*.c
        ${CMAKE_SOURCE_DIR}/src/gromacs/linearalgebra/gmx_blas/*.c
        ${CMAKE_SOURCE_DIR}/src/gromacs/linearalgebra/gmx_lapack/*.c
        ${CMAKE_SOURCE_DIR}/src/contrib/*.c
        ${CMAKE_SOURCE_DIR}/src/contrib/*.cpp
        ${CMAKE_SOURCE_DIR}/src/contrib/*.cu
        ${CMAKE_SOURCE_DIR}/src/external/*.c
        ${CMAKE_SOURCE_DIR}/src/external/*.cpp
        ${CMAKE_SOURCE_DIR}/src/external/*.cu
        ${CMAKE_SOURCE_DIR}/src/gromacs/selection/scanner.cpp
        ${CMAKE_SOURCE_DIR}/src/gromacs/selection/parser.cpp
        ${CMAKE_SOURCE_DIR}/src/gromacs/gpu_utils/ocl_compiler.cpp
        )
    list(REMOVE_ITEM _inputfiles ${_files_to_ignore})

    # Set flags for cppcheck
    set(_outputext txt)
    set(_outputopt "--template={file}:{line}:{id}: {severity}: {message}")
    if (CPPCHECK_XML_OUTPUT)
        set(_outputext xml)
        set(_outputopt --xml --xml-version=2)
    endif()
    set(_common_flags
        --enable=style -DLINUX -DHAVE_UNISTD_H
        -I src
        -I src/external/thread_mpi/include
        -I src/external/tng_io/include
        -I ${CMAKE_BINARY_DIR}/src
        --quiet
        --inline-suppr
        ${_outputopt})
    set(_c_flags
        --suppress=variableScope
        --suppress=unnecessaryForwardDeclaration
        --suppress=unusedVariable
        --suppress=unreadVariable
        --suppress=unusedStructMember
        --suppress=invalidscanf
        --suppress=sizeofCalculation
        --suppress=invalidscanf_libc
        --suppress=*:src/external/Random123-1.08/include/Random123/features/compilerfeatures.h
        ) 
    set(_cxx_flags
        --language=c++ # particularly useful for CUDA files
        -D__cplusplus
        --suppress=variableScope
        --suppress=unnecessaryForwardDeclaration
        --suppress=memsetClassFloat  #we assume IEEE754
        --suppress=invalidscanf_libc #seems only important for security on non-std libc
        --suppress=invalidscanf      #same as last (style and portability checker have the same warning)
        --suppress=passedByValue:src/gromacs/simd/tests/*
        --suppress=redundantAssignment:src/gromacs/simd/simd_math.h #seems to be a bug in cppcheck
        --suppress=noExplicitConstructor # can't be selective about this, see http://sourceforge.net/p/cppcheck/discussion/general/thread/db1e4ba7/
        --suppress=unusedStructMember:src/gromacs/onlinehelp/tests/helpmanager.cpp
        --suppress=unusedStructMember:src/gromacs/commandline/cmdlinehelpmodule.cpp
        --suppress=unusedStructMember:src/gromacs/selection/selhelp.cpp
        --suppress=redundantPointerOp:src/gromacs/fileio/gmxfio-xdr.cpp
        --suppress=passedByValue # See comment below
        --suppress=shiftTooManyBits:src/gromacs/gpu_utils/gpu_utils.cu # CUDA kernel launch false positive
        )
        # Passing non-trivial objects by value is rarely a problem for
        # GROMACS in performance-sensitive code, and shouldn't be
        # enforced for types that are intended to be used like value
        # types (e.g. SIMD wrapper types, ArrayRef) , nor for
        # move-enabled types. cppcheck isn't sensitive to these
        # subtleties yet.

    # This list will hold the list of all files with cppcheck errors
    # (one per input file)
    set(_filelist)

    foreach (_filename ${_inputfiles})
        set(_target_name cppcheck-${_filename}.${_outputext})
        string(REPLACE "/" "_" _target_name ${_target_name})
        list(APPEND _filelist ${_target_name})
        if (_filename MATCHES "\\.cpp$" OR _filename MATCHES "\\.cu$")
            set(_lang CXX)
            set(_lang_flags ${_cxx_flags})
        else()
            set(_lang C)
            set(_lang_flags ${_c_flags})
        endif()
        add_custom_command(
            OUTPUT ${_target_name}
            COMMAND ${CPPCHECK_EXECUTABLE} ${_common_flags} ${_lang_flags}
                ${_filename} 2> ${CMAKE_CURRENT_BINARY_DIR}/${_target_name}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            IMPLICIT_DEPENDS ${_lang} ${CMAKE_SOURCE_DIR}/${_filename}
            COMMENT "Running cppcheck for ${_filename}"
            VERBATIM)
        if (NOT CPPCHECK_XML_OUTPUT)
            add_custom_command(
                OUTPUT ${_target_name}
                COMMAND cat ${CMAKE_CURRENT_BINARY_DIR}/${_target_name}
                VERBATIM APPEND)
        endif()
    endforeach()
    if (NOT CPPCHECK_XML_OUTPUT)
        set(_target_name cppcheck-errors.${_outputext})
        add_custom_command(
            OUTPUT ${_target_name}
            COMMAND cat ${_filelist} > ${_target_name}
            DEPENDS ${_filelist}
            COMMENT "Combining cppcheck results"
            VERBATIM)
        list(APPEND _filelist ${_target_name})
    endif()
    add_custom_target(cppcheck DEPENDS ${_filelist})
else()
    add_custom_target(cppcheck
        COMMAND ${CMAKE_COMMAND} -E echo
            "cppcheck was not found by CMake. Rerun CMake specifying CPPCHECK_EXECUTABLE."
        VERBATIM)
endif()
