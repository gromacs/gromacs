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

if (NOT GMX_BUILD_UNITTESTS)
    gmx_add_missing_tests_notice("Unit tests have not been run. You need to set GMX_BUILD_UNITTESTS=ON if you want to build and run them.")
    return()
endif()

include_directories(BEFORE SYSTEM ${GMOCK_INCLUDE_DIRS})
set(TESTUTILS_SOURCES
    cmdlinetest.cpp
    conftest.cpp
    filematchers.cpp
    interactivetest.cpp
    loggertest.cpp
    mpi-printer.cpp
    mpitest.cpp
    refdata.cpp
    refdata-xml.cpp
    stdiohelper.cpp
    stringtest.cpp
    testasserts.cpp
    testfilemanager.cpp
    testfileredirector.cpp
    testinit.cpp
    testoptions.cpp
    TextSplitter.cpp
    textblockmatchers.cpp
    xvgtest.cpp
    )

if(NOT HAVE_TINYXML2)
    list(APPEND TESTUTILS_SOURCES ../external/tinyxml2/tinyxml2.cpp)
endif()

add_library(testutils STATIC ${UNITTEST_TARGET_OPTIONS} ${TESTUTILS_SOURCES})
set(TESTUTILS_LIBS testutils)
set_property(TARGET testutils APPEND PROPERTY COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}")
set_property(TARGET testutils APPEND PROPERTY COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
target_link_libraries(testutils libgromacs ${GMX_COMMON_LIBRARIES} ${GMOCK_LIBRARIES})

if(HAVE_TINYXML2)
    include_directories(SYSTEM ${TinyXML2_INCLUDE_DIR})
    target_link_libraries(testutils ${TinyXML2_LIBRARIES})
else()
    include_directories(BEFORE SYSTEM "../external/tinyxml2")
endif()

# TODO Use gmx_add_missing_tests_notice() instead of the messages below.
set(GMX_CAN_RUN_MPI_TESTS 1)
if (GMX_MPI)
    set(_an_mpi_variable_had_content 0)
    foreach(VARNAME MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS)
        # These variables need a valid value for the test to run
        # and pass, but conceivably any of them might be valid
        # with arbitrary (including empty) content. They can't be
        # valid if they've been populated with the CMake
        # find_package magic suffix/value "NOTFOUND", though.
        if (${VARNAME} MATCHES ".*NOTFOUND")
            message(STATUS "CMake variable ${VARNAME} was not detected to be a valid value. To test GROMACS correctly, check the advice in the install guide.")
            set(GMX_CAN_RUN_MPI_TESTS 0)
        endif()
        if (NOT VARNAME STREQUAL MPIEXEC AND ${VARNAME})
            set(_an_mpi_variable_had_content 1)
        endif()
    endforeach()
    if(_an_mpi_variable_had_content AND NOT MPIEXEC)
        message(STATUS "CMake variable MPIEXEC must have a valid value if one of the other related MPIEXEC variables does. To test GROMACS correctly, check the advice in the install guide.")
        set(GMX_CAN_RUN_MPI_TESTS 0)
    endif()
elseif (NOT GMX_THREAD_MPI)
    set(GMX_CAN_RUN_MPI_TESTS 0)
endif()

set(TESTUTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TESTUTILS_DIR ${TESTUTILS_DIR} PARENT_SCOPE)
set(TESTUTILS_LIBS ${TESTUTILS_LIBS} PARENT_SCOPE)
set(GMX_CAN_RUN_MPI_TESTS ${GMX_CAN_RUN_MPI_TESTS} PARENT_SCOPE)

add_subdirectory(tests)
