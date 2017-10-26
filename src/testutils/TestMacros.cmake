#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

include(CMakeParseArguments)

function (gmx_add_unit_test_object_library NAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        include_directories(BEFORE SYSTEM ${GMOCK_INCLUDE_DIRS})
        add_library(${NAME} OBJECT ${UNITTEST_TARGET_OPTIONS} ${ARGN})
        set_property(TARGET ${NAME} APPEND PROPERTY COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}")
        set_property(TARGET ${NAME} APPEND PROPERTY COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
    endif()
endfunction ()

function (gmx_add_gtest_executable EXENAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        set(_options MPI HARDWARE_DETECTION)
        cmake_parse_arguments(ARG "${_options}" "" "" ${ARGN})
        set(_source_files ${ARG_UNPARSED_ARGUMENTS})

        file(RELATIVE_PATH _input_files_path ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
        set(_temporary_files_path "${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary")
        file(MAKE_DIRECTORY ${_temporary_files_path})
        # Note that the quotation marks in the next line form part of
        # the defined symbol, so that the macro replacement in the
        # source file is as a string.
        # These are only needed for unittest_main.cpp, but for simplicity used
        # for the whole target (since there may be multiple executables in the
        # same directory, it is not straightforward to use a source file
        # property).
        set(EXTRA_COMPILE_DEFINITIONS
            TEST_DATA_PATH="${_input_files_path}"
            TEST_TEMP_PATH="${_temporary_files_path}")
        if (ARG_MPI)
            list(APPEND EXTRA_COMPILE_DEFINITIONS
                 TEST_USES_MPI=true)
        endif()
        if (ARG_HARDWARE_DETECTION)
            list(APPEND EXTRA_COMPILE_DEFINITIONS
                 TEST_USES_HARDWARE_DETECTION=true)
        endif()

        include_directories(BEFORE SYSTEM ${GMOCK_INCLUDE_DIRS})
        add_executable(${EXENAME} ${UNITTEST_TARGET_OPTIONS}
            ${_source_files} ${TESTUTILS_DIR}/unittest_main.cpp)
        target_link_libraries(${EXENAME}
            ${TESTUTILS_LIBS} libgromacs ${GMOCK_LIBRARIES}
            ${GMX_COMMON_LIBRARIES} ${GMX_EXE_LINKER_FLAGS} ${GMX_STDLIB_LIBRARIES})
        set_property(TARGET ${EXENAME}
            APPEND PROPERTY COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
        set_property(TARGET ${EXENAME}
            APPEND PROPERTY COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}")
        set_property(TARGET ${EXENAME}
            APPEND PROPERTY COMPILE_DEFINITIONS "${EXTRA_COMPILE_DEFINITIONS}")
    endif()
endfunction()

# Use this function with MPI_RANKS <N> INTEGRATION_TEST to register a test
# binary as an integration test that requires MPI. The intended number of MPI
# ranks is also passed
#
# TODO When a test case needs it, generalize the MPI_RANKS mechanism so
# that ctest can run the test binary over a range of numbers of MPI
# ranks.
function (gmx_register_gtest_test NAME EXENAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        set(_options INTEGRATION_TEST)
        set(_one_value_args MPI_RANKS)
        cmake_parse_arguments(ARG "${_options}" "${_one_value_args}" "" ${ARGN})
        set(_xml_path ${CMAKE_BINARY_DIR}/Testing/Temporary/${NAME}.xml)
        set(_labels GTest)
        set(_timeout 30)
        if (ARG_INTEGRATION_TEST)
            list(APPEND _labels IntegrationTest)
            set(_timeout 120)
            gmx_get_test_prefix_cmd(_prefix_cmd IGNORE_LEAKS)
        else()
            list(APPEND _labels UnitTest)
            gmx_get_test_prefix_cmd(_prefix_cmd)
        endif()
        set(_cmd ${_prefix_cmd} $<TARGET_FILE:${EXENAME}>)
        if (ARG_MPI_RANKS)
            if (NOT GMX_CAN_RUN_MPI_TESTS)
                return()
            endif()
            list(APPEND _labels MpiTest)
            if (GMX_MPI)
                set(_cmd
                    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ARG_MPI_RANKS}
                    ${MPIEXEC_PREFLAGS} ${_cmd} ${MPIEXEC_POSTFLAGS})
            elseif (GMX_THREAD_MPI)
                list(APPEND _cmd -ntmpi ${ARG_MPI_RANKS})
            endif()
        endif()
        add_test(NAME ${NAME}
                 COMMAND ${_cmd} --gtest_output=xml:${_xml_path})
        set_tests_properties(${NAME} PROPERTIES LABELS "${_labels}")
        set_tests_properties(${NAME} PROPERTIES TIMEOUT ${_timeout})
        add_dependencies(tests ${EXENAME})
    endif()
endfunction ()

function (gmx_add_unit_test NAME EXENAME)
    gmx_add_gtest_executable(${EXENAME} ${ARGN})
    gmx_register_gtest_test(${NAME} ${EXENAME})
endfunction()

function (gmx_add_mpi_unit_test NAME EXENAME RANKS)
    if (GMX_MPI OR (GMX_THREAD_MPI AND GTEST_IS_THREADSAFE))
        gmx_add_gtest_executable(${EXENAME} MPI ${ARGN})
        gmx_register_gtest_test(${NAME} ${EXENAME} MPI_RANKS ${RANKS})
    endif()
endfunction()
