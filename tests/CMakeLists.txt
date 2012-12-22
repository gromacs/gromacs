#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
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
enable_testing()
add_test(TestExec_mdrun-h ../src/kernel/mdrun${GMX_BINARY_SUFFIX} -h)
if(NOT EXISTS ${REGRESSIONTEST_PATH})
    if(REGRESSIONTEST_DOWNLOAD)
        if(GITCOMMAND)
            #https is most likely to work. Otherwise the user should manually download
            #git download can be replaced with TGZ when TGZ is available.
            execute_process(COMMAND ${GITCOMMAND} 
                clone https://gerrit.gromacs.org/regressiontests
                RESULT_VARIABLE DOWNLOAD_FAILED)
            if(DOWNLOAD_FAILED)
                message(WARNING "Download of tests failed. Download manually and set REGRESSIONTEST_PATH if you want to run the tests.")
            else()
                set(REGRESSIONTEST_PATH ${CMAKE_BINARY_DIR}/regressiontests CACHE 
                    STRING "Path to regressiontests" FORCE)
                mark_as_advanced(REGRESSIONTEST_PATH)
            endif()
        else()
            message(WARNING "Please set GITCOMMAND or REGRESSIONTEST_PATH")
        endif()
    else()
        message(WARNING "To run the regressiontests, you need to either set REGRESSIONTEST_PATH to a valid path, or set REGRESSIONTEST_DOWNLOAD to yes. The later will download them automatically if you have git available.")
    endif()
endif()
set(ARGS "")
if(GMX_DOUBLE)
    set(ARGS ${ARGS} -double)
endif()
if(GMX_LIB_MPI AND NOT MPIEXEC) #autodetection failed or CC=mpicc was used
    message(WARNING "Please set MPIEXEC. Otherwise mpirun is assumed for runnings tests.")
endif()
if(GMX_LIB_MPI)
    set(GMX_TEST_NUMBER_PROCS 8 CACHE STRING "Number of processors used for testing")
    mark_as_advanced(GMX_TEST_NUMBER_PROCS)
    set(ARGS ${ARGS} -np ${GMX_TEST_NUMBER_PROCS})
    if(MPIEXEC)
        set(ARGS ${ARGS} -mpirun ${MPIEXEC})
    endif()
endif()
if(GMX_BLUEGENE)
    set(ARGS ${ARGS} -bluegene)
endif()
if(GMX_BINARY_SUFFIX)
    set(ARGS -suffix ${GMX_BINARY_SUFFIX})
endif()
#crosscompile is only used to disable checking whether binaries work
#given that we know they are there and that mpirun might not be exectuable (e.g. Cray)
#we enable it.
set(ARGS ${ARGS} -crosscompile)

if(GMX_NATIVE_WINDOWS)
    set(DIR_SEPERATOR "\\;")
    #replacing \ with / shouldn't be neccessary. But otherwise "..\bin\;c:\.." gets turned into "...\bin\\c:\.."
    #don't know why and don't have a better workaround. This workaround doesn't hurt.
    string(REPLACE "\\" "/" PATH "$ENV{PATH}")
    #protect ; (don't treat as list)
    string(REPLACE ";" "\\;" PATH "${PATH}")
else()
    set(DIR_SEPERATOR ":")
    set(PATH "$ENV{PATH}")
endif()

foreach(FOLDER kernel tools gmxlib mdlib) #lib folders might be needed for e.g. DLLs
    #for GMX paths native ("\") is needed so that GMXLIB detection works
    file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/src/${FOLDER}" DIR)
    set(PATH "${DIR}${DIR_SEPERATOR}${PATH}")
endforeach()
#currently not testing tools because they don't contain any useful tests
foreach(subtest simple complex kernel freeenergy pdb2gmx)
    add_test(NAME GmxTest-${subtest} 
        WORKING_DIRECTORY "${REGRESSIONTEST_PATH}"
        #windows requires the command to be perl and not the script
        COMMAND perl "${REGRESSIONTEST_PATH}/gmxtest.pl" ${subtest} -noverbose 
            -nosuffix ${ARGS})
    set_tests_properties(GmxTest-${subtest} PROPERTIES ENVIRONMENT "PATH=${PATH}"
        REQUIRED_FILES "${REGRESSIONTEST_PATH}/gmxtest.pl")
endforeach()
