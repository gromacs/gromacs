#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

function (gmx_add_unit_test_object_library NAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        include_directories(BEFORE SYSTEM ${GMOCK_INCLUDE_DIRS})
        add_library(${NAME} OBJECT ${UNITTEST_TARGET_OPTIONS} ${ARGN})
        set_property(TARGET ${NAME} APPEND PROPERTY COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}")
        set_property(TARGET ${NAME} APPEND PROPERTY COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
    endif()
endfunction ()

function (gmx_build_unit_test NAME EXENAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        include_directories(BEFORE SYSTEM ${GMOCK_INCLUDE_DIRS})
        add_executable(${EXENAME} ${UNITTEST_TARGET_OPTIONS} ${ARGN} ${TESTUTILS_DIR}/unittest_main.cpp)
        set_property(TARGET ${EXENAME} APPEND PROPERTY COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}")
        set_property(TARGET ${EXENAME} APPEND PROPERTY COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
        target_link_libraries(${EXENAME} ${TESTUTILS_LIBS} libgromacs ${GMOCK_LIBRARIES} ${GMX_EXE_LINKER_FLAGS})
        file(RELATIVE_PATH _input_files_path ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
        set(_temporary_files_path "${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary")
        file(MAKE_DIRECTORY ${_temporary_files_path})
        # Note that the quotation marks in the next line form part of
        # the defined symbol, so that the macro replacement in the
        # source file is as a string.
        set(EXTRA_COMPILE_DEFINITIONS
            TEST_DATA_PATH="${_input_files_path}"
            TEST_TEMP_PATH="${_temporary_files_path}")

        set_property(TARGET ${EXENAME} APPEND PROPERTY COMPILE_DEFINITIONS "${EXTRA_COMPILE_DEFINITIONS}")
    endif()
endfunction ()

function (gmx_register_unit_test NAME EXENAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        add_test(NAME ${NAME}
                 COMMAND ${EXENAME} --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/Temporary/${NAME}.xml)
        set_tests_properties(${NAME} PROPERTIES LABELS "GTest;UnitTest")
        add_dependencies(tests ${EXENAME})
    endif()
endfunction ()

# Use this function to register a test binary as an integration test
function (gmx_register_integration_test NAME EXENAME)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        add_test(NAME ${NAME}
                 COMMAND ${EXENAME} --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/Temporary/${NAME}.xml)
        set_tests_properties(${testname} PROPERTIES LABELS "IntegrationTest")
        add_dependencies(tests ${EXENAME})

        # GMX_EXTRA_LIBRARIES might be needed for mdrun integration tests at
        # some point.
        # target_link_libraries(${EXENAME} ${GMX_EXTRA_LIBRARIES})
    endif()
endfunction ()

# Use this function to register a test binary as an integration test
# that requires MPI. The intended number of MPI ranks is also passed
#
# TODO When a test case needs it, generalize the NUMPROC mechanism so
# that ctest can run the test binary over a range of numbers of MPI
# ranks.
function (gmx_register_mpi_integration_test NAME EXENAME NUMPROC)
    if (GMX_BUILD_UNITTESTS AND BUILD_TESTING)
        if (GMX_MPI)
            foreach(VARNAME MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS)
                # These variables need a valid value for the test to run
                # and pass, but conceivably any of them might be valid
                # with arbitrary (including empty) content. They can't be
                # valid if they've been populated with the CMake
                # find_package magic suffix/value "NOTFOUND", though.
                if (${VARNAME} MATCHES ".*NOTFOUND")
                    message(STATUS "CMake variable ${VARNAME} was not detected to be a valid value. To test GROMACS correctly, check the advice in the install guide.")
                    set(_cannot_run_mpi_tests 1)
                endif()
                if (NOT VARNAME STREQUAL MPIEXEC AND ${VARNAME})
                    set(_an_mpi_variable_had_content 1)
                endif()
            endforeach()
            if(_an_mpi_variable_had_content AND NOT MPIEXEC)
                message(STATUS "CMake variable MPIEXEC must have a valid value if one of the other related MPIEXEC variables does. To test GROMACS correctly, check the advice in the install guide.")
                set(_cannot_run_mpi_tests 1)
            endif()
            if(NOT _cannot_run_mpi_tests)
                add_test(NAME ${NAME}
                    COMMAND
                    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC}
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${EXENAME}> ${MPIEXEC_POSTFLAGS}
                    --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/Temporary/${NAME}.xml
                    )
                set_tests_properties(${testname} PROPERTIES LABELS "MpiIntegrationTest")
                add_dependencies(tests ${EXENAME})
            endif()

            # GMX_EXTRA_LIBRARIES might be needed for mdrun integration tests at
            # some point.
            # target_link_libraries(${EXENAME} ${GMX_EXTRA_LIBRARIES})
        elseif(GMX_THREAD_MPI)
            add_test(NAME ${NAME}
                COMMAND
                $<TARGET_FILE:${EXENAME}> -nt ${NUMPROC}
                --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/Temporary/${NAME}.xml
                )
            set_tests_properties(${testname} PROPERTIES LABELS "MpiIntegrationTest")
            add_dependencies(tests ${EXENAME})

            # GMX_EXTRA_LIBRARIES might be needed for mdrun integration tests at
            # some point.
            # target_link_libraries(${EXENAME} ${GMX_EXTRA_LIBRARIES})
        endif()
    endif()
endfunction ()

function (gmx_add_unit_test NAME EXENAME)
    gmx_build_unit_test(${NAME} ${EXENAME} ${ARGN})
    gmx_register_unit_test(${NAME} ${EXENAME})
endfunction()
