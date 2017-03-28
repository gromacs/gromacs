#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

set(REGRESSIONTEST_PATH "" CACHE PATH "Directory containing regressiontests")
mark_as_advanced(REGRESSIONTEST_PATH)
option(REGRESSIONTEST_DOWNLOAD
    "Automatically download regressiontests. Tests can be run with ctest." OFF)

if(REGRESSIONTEST_DOWNLOAD)
    if (NOT SOURCE_IS_SOURCE_DISTRIBUTION)
        set(REGRESSIONTEST_URL http://gerrit.gromacs.org/snapshot/${REGRESSIONTEST_BRANCH})
        # REGRESSIONTEST_PATH for dev trees is set later based on the dirname found in the tar
    else()
        set(REGRESSIONTEST_URL http://gerrit.gromacs.org/download/regressiontests-${REGRESSIONTEST_VERSION}.tar.gz)
        set(REGRESSIONTEST_PATH
            "${CMAKE_CURRENT_BINARY_DIR}/regressiontests-${REGRESSIONTEST_VERSION}"
            CACHE PATH "Path to auto-downloaded regressiontests" FORCE)
    endif()
    set(REGRESSIONTEST_FILE "${CMAKE_CURRENT_BINARY_DIR}/regressiontests.tgz")
    message("Downloading: ${REGRESSIONTEST_URL}")
    file(DOWNLOAD ${REGRESSIONTEST_URL} "${REGRESSIONTEST_FILE}" SHOW_PROGRESS STATUS status LOG log)
    list(GET status 0 status_code)
    list(GET status 1 status_string)

    if(NOT status_code EQUAL 0)
        message(FATAL_ERROR "error: downloading '${REGRESSIONTEST_URL}' failed
status_code: ${status_code}
status_string: ${status_string}
log: ${log}")
    endif()
    if (SOURCE_IS_SOURCE_DISTRIBUTION)
        file(MD5 ${REGRESSIONTEST_FILE} COMPUTED_MD5SUM)
        if(NOT ${REGRESSIONTEST_MD5SUM} STREQUAL ${COMPUTED_MD5SUM})
            message(FATAL_ERROR "Download of regressiontests failed. Expected MD5 of ${REGRESSIONTEST_MD5SUM} but download has ${COMPUTED_MD5SUM}")
        endif()
    else()
        # Extract the REGRESSIONTEST_PATH from the tar when using dev tree.
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar tf "${REGRESSIONTEST_FILE}"
                RESULT_VARIABLE _tar_tf_res
                OUTPUT_VARIABLE _tar_tf_out)
        if (${_tar_tf_res} EQUAL 0)
            string(REGEX REPLACE "/\n.*$" "" _regressiontest_dir "${_tar_tf_out}")
            set(REGRESSIONTEST_PATH "${CMAKE_CURRENT_BINARY_DIR}/${_regressiontest_dir}"
                CACHE PATH "Path to auto-downloaded regressiontests" FORCE)
        else()
            message(FATAL_ERROR "Failed to list the contents of the downloaded tarball: ${REGRESSIONTEST_FILE}")
        endif()
    endif()

    file(REMOVE_RECURSE "${REGRESSIONTEST_PATH}") #delete potential prior folder
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf "${REGRESSIONTEST_FILE}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    if(NOT EXISTS ${REGRESSIONTEST_PATH}/gmxtest.pl)
        message(FATAL_ERROR "Download incorrect. Doesn't contain required gmxtest.pl")
    endif()
    set(REGRESSIONTEST_DOWNLOAD OFF CACHE BOOL "Tests already downloaded. Set to yes to download again" FORCE)
endif()

if(REGRESSIONTEST_PATH AND (CMAKE_CROSSCOMPILING OR CMAKE_CONFIGURATION_TYPES OR GMX_BUILD_MDRUN_ONLY))
    # TODO: It would be nicer to do these checks before potentially downloading the tests.
    # Cross-compiling toolchains require us to compile both front-end and
    # back-end binaries to run gmxtest.pl.
    # Testing an mdrun-only builds require supporting binaries from a full build
    message(WARNING
        "With cross-compiling, multi-configuration generators (e.g. Visual Studio), or with mdrun-only builds, running regressiontests from build system is not supported. Please run gmxtest.pl directly.")
    set(REGRESSIONTEST_PATH OFF CACHE BOOL
        "With cross-compiling or multi-configuration generators, running regressiontests from build system is not supported." FORCE)
endif()

if(REGRESSIONTEST_PATH)
    if(NOT EXISTS ${REGRESSIONTEST_PATH}/gmxtest.pl)
        message(FATAL_ERROR
            "REGRESSIONTEST_PATH invalid. The path needs to contain gmxtest.pl.")
    endif()

    # gmxtests target builds all binaries required for running gmxtest
    add_custom_target(gmxtests DEPENDS gmx)
    add_dependencies(run-ctest gmxtests)
    add_dependencies(run-ctest-nophys gmxtests)

    if(GMX_DOUBLE)
        list(APPEND ARGS -double)
    endif()
    if(GMX_LIB_MPI AND NOT MPIEXEC) # find_program failed
        message(WARNING
            "Please set MPIEXEC. Otherwise mpirun is assumed for running tests.")
    endif()
    if(GMX_LIB_MPI)
        set(GMX_TEST_NUMBER_PROCS 8 CACHE STRING "Number of processors used for testing")
        mark_as_advanced(GMX_TEST_NUMBER_PROCS)
        list(APPEND ARGS -np ${GMX_TEST_NUMBER_PROCS})
        if(MPIEXEC)
            list(APPEND ARGS -mpirun ${MPIEXEC})
        endif()
        #We should use MPIEXEC_NUMPROC_FLAG but gmxtest.pl doesn't let us pass it
    endif()
    if(GMX_BINARY_SUFFIX)
        list(APPEND ARGS -suffix ${GMX_BINARY_SUFFIX})
    endif()
    #crosscompile is only used to disable checking whether binaries work
    #given that we know they are there and that mdrun might not be exectuable
    #(e.g. Cray) we enable it.
    list(APPEND ARGS -crosscompile)

    set(REGRESSIONTEST_EXTRA_ARGS "" CACHE STRING 
        "Extra arguments passed to gmxtest")
    mark_as_advanced(REGRESSIONTEST_EXTRA_ARGS)
    list(APPEND ARGS ${REGRESSIONTEST_EXTRA_ARGS})

    list(APPEND ARGS -noverbose -nosuffix)

    if(GMX_NATIVE_WINDOWS)
        set(PATH_SEPARATOR "\\;")
        #replacing \ with / shouldn't be neccessary. But otherwise "..\bin\;c:\.."
        #gets turned into "...\bin\\c:\.." don't know why and don't have a better
        #workaround. This workaround doesn't hurt.
        string(REPLACE "\\" "/" PATH "$ENV{PATH}")
        #protect ; (don't treat as list)
        string(REPLACE ";" "\\;" PATH "${PATH}")
    else()
        set(PATH_SEPARATOR ":")
        set(PATH "$ENV{PATH}")
    endif()

    foreach(FOLDER bin lib) #lib folders might be needed for
        #e.g. DLLs. For GMX paths native ("\") is needed for GMXLIB detection
        file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/${FOLDER}" DIR)
        set(PATH "${DIR}${PATH_SEPARATOR}${PATH}")
    endforeach()

    find_program(PERL_EXECUTABLE NAMES "perl")
    mark_as_advanced(PERL_EXECUTABLE)

    if (NOT PERL_EXECUTABLE)
        message(FATAL_ERROR "Perl not found. Install perl, set PERL_EXECUTABLE to the perl location, or unset REGRESSIONTEST_PATH to disable testing.")
    endif()

    #currently not testing tools because they don't contain any useful tests
    foreach(subtest simple complex kernel freeenergy pdb2gmx rotation)
        add_test(NAME regressiontests/${subtest}
            #windows requires the command to be perl and not the script
            COMMAND perl "${REGRESSIONTEST_PATH}/gmxtest.pl" ${subtest} ${ARGS})
        set_tests_properties(regressiontests/${subtest} PROPERTIES
            ENVIRONMENT "PATH=${PATH}")
    endforeach()
else()
    gmx_add_missing_tests_notice_nophys("Regression tests have not been run. If you want to run them from the build system, get the correct version of the regression tests package and set REGRESSIONTEST_PATH in CMake to point to it, or set REGRESSIONTEST_DOWNLOAD=ON.")
endif()


#
# Create physical validation test in analogy to regression tests
#

set(PHYSVALTEST_PATH "" CACHE PATH "Directory containing physical validation tests")
mark_as_advanced(PHYSVALTEST_PATH)

# TODO: Discuss distribution of physical validation suite. Could be either with the source, as
#       done currently (in admin/builds/physicalvalidation), or via download, as for the
#       regression tests. As it is distributed with the source currently, let's copy the
#       directory and set the path explicitly. The copying is necessary to make sure that to
#       builds from the same source do not interfere with each others while testing.
#       The path can still be changed in the cache, if the need arises.
file(COPY ${CMAKE_SOURCE_DIR}/admin/builds/physicalvalidation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
set(PHYSVALTEST_PATH "${CMAKE_CURRENT_BINARY_DIR}/physicalvalidation" CACHE PATH "Directory containing physical validation tests" FORCE)

# TODO: Ensure that the following is also true for physical validation tests (copied from regtests)!
#
if(PHYSVALTEST_PATH AND (CMAKE_CROSSCOMPILING OR CMAKE_CONFIGURATION_TYPES OR GMX_BUILD_MDRUN_ONLY))
    # TODO: It would be nicer to do these checks before potentially downloading the tests.
    # Cross-compiling toolchains require us to compile both front-end and
    # back-end binaries to run gmxtest.pl.
    # Testing an mdrun-only builds require supporting binaries from a full build
    message(WARNING
        "With cross-compiling, multi-configuration generators (e.g. Visual Studio), or with mdrun-only builds, running physical validation tests from build system is not supported. Please run physicalvalidation.py directly.")
    set(PHYSVALTEST_PATH OFF CACHE BOOL
        "With cross-compiling or multi-configuration generators, running physical validation tests from build system is not supported." FORCE)
endif()

if(PHYSVALTEST_PATH)
    #
    # Check if path valid
    #
    if(NOT EXISTS ${PHYSVALTEST_PATH}/gmx_physicalvalidation.py)
        message(FATAL_ERROR
            "PHYSVALTEST_PATH invalid. The path needs to contain gmx_physicalvalidation.py.")
    endif()

    #
    # Build target
    #
    # physvaltests target builds binaries required for running physical validation
    add_custom_target(physvaltests DEPENDS gmx)
    add_dependencies(run-ctest-phys physvaltests)
    add_dependencies(run-ctest physvaltests)

    #
    # Determine arguments passed to physicalvalidation.py
    #
    if(GMX_DOUBLE)
        list(APPEND PARGS --double)
    endif()

    # TODO: Handle MPI in physicalvalidation.py
    #
    # if(GMX_LIB_MPI AND NOT MPIEXEC) # find_program failed
    #     message(WARNING
    #         "Please set MPIEXEC. Otherwise mpirun is assumed for running tests.")
    # endif()
    # if(GMX_LIB_MPI)
    #     set(GMX_TEST_NUMBER_PROCS 8 CACHE STRING "Number of processors used for testing")
    #     mark_as_advanced(GMX_TEST_NUMBER_PROCS)
    #     list(APPEND ARGS -np ${GMX_TEST_NUMBER_PROCS})
    #     if(MPIEXEC)
    #         list(APPEND ARGS -mpirun ${MPIEXEC})
    #     endif()
    #     #We should use MPIEXEC_NUMPROC_FLAG but gmxtest.pl doesn't let us pass it
    # endif()
    if(GMX_LIB_MPI)
        message(FATAL_ERROR
            "Physical validation using MPI not supported.")
    endif()

    if(GMX_BINARY_SUFFIX)
        list(APPEND ARGS --suffix ${GMX_BINARY_SUFFIX})
    endif()
    
    # TODO: Understand if this is something we need
    #
    # #crosscompile is only used to disable checking whether binaries work
    # #given that we know they are there and that mdrun might not be exectuable
    # #(e.g. Cray) we enable it.
    # list(APPEND ARGS -crosscompile)

    set(PHYSVALTEST_EXTRA_ARGS "" CACHE STRING 
        "Extra arguments passed to phystest")
    mark_as_advanced(PHYSVALTEST_EXTRA_ARGS)
    list(APPEND PARGS ${PHYSVALTEST_EXTRA_ARGS})

    # add extra flags
    #
    # list(APPEND PARGS --)

    #
    # The following lines are directly copied from regression tests.
    # They seem to be applicable to physical validation tests as well.
    #
    if(GMX_NATIVE_WINDOWS)
        set(PATH_SEPARATOR "\\;")
        #replacing \ with / shouldn't be neccessary. But otherwise "..\bin\;c:\.."
        #gets turned into "...\bin\\c:\.." don't know why and don't have a better
        #workaround. This workaround doesn't hurt.
        string(REPLACE "\\" "/" PATH "$ENV{PATH}")
        #protect ; (don't treat as list)
        string(REPLACE ";" "\\;" PATH "${PATH}")
    else()
        set(PATH_SEPARATOR ":")
        set(PATH "$ENV{PATH}")
    endif()

    foreach(FOLDER bin lib) #lib folders might be needed for
                            #e.g. DLLs. For GMX paths native ("\") is needed for GMXLIB detection
        file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/${FOLDER}" DIR)
        set(PATH "${DIR}${PATH_SEPARATOR}${PATH}")
    endforeach()
    #
    # End copied from regression tests.
    #

    #
    # Test if we find python3
    #
    find_program(PYTHON_EXECUTABLE NAMES "python3")
    mark_as_advanced(PYTHON_EXECUTABLE)
    if (NOT PYTHON_EXECUTABLE)
        message(FATAL_ERROR "Python3 not found. Install python3, set PYTHON_EXECUTABLE to the python3 location, or unset PHYSVALTEST_PATH to disable testing.")
    endif()

    #
    # Hook in our own tests
    #
    foreach(subtest pass integratorconvergence)
        add_test(NAME physvaltests/${subtest}
            # For reg tests: "windows requires the command to be perl and not the script"
            COMMAND python3 "${PHYSVALTEST_PATH}/gmx_physicalvalidation.py" ${subtest} ${PARGS})
        set_tests_properties(physvaltests/${subtest} PROPERTIES
            ENVIRONMENT "PATH=${PATH}")
    endforeach()
else()
    gmx_add_missing_tests_notice_phys("Physical validation tests have not been run. If you want to run them from the build system, get the correct version of the physical validation tests package and set PHYSVALTEST_PATH in CMake to point to it.") #, or set PHYSVALTEST_DOWNLOAD=ON.")
endif()


gmx_create_missing_tests_notice_target()
include(CppCheck.cmake)
