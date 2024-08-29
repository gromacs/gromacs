#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

set(REGRESSIONTEST_PATH "" CACHE PATH "Directory containing regressiontests")
mark_as_advanced(REGRESSIONTEST_PATH)
option(REGRESSIONTEST_DOWNLOAD
    "Automatically download regressiontests. Tests can be run with ctest." OFF)

if(REGRESSIONTEST_DOWNLOAD)
    if (NOT SOURCE_IS_SOURCE_DISTRIBUTION)
        set(REGRESSIONTEST_URL https://gitlab.com/gromacs/gromacs-regressiontests/-/archive/${REGRESSIONTEST_BRANCH}/gromacs-regressiontests-${REGRESSIONTEST_BRANCH}.tar.gz)
        # REGRESSIONTEST_PATH for dev trees is set later based on the dirname found in the tar
    else()
        set(REGRESSIONTEST_URL https://ftp.gromacs.org/regressiontests/regressiontests-${REGRESSIONTEST_VERSION}.tar.gz)
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

if(REGRESSIONTEST_PATH AND (CMAKE_CROSSCOMPILING OR CMAKE_CONFIGURATION_TYPES))
    # TODO: It would be nicer to do these checks before potentially downloading the tests.
    # Cross-compiling toolchains require us to compile both front-end and
    # back-end binaries to run gmxtest.pl.
    message(WARNING
        "With cross-compiling or multi-configuration generators (e.g. Visual Studio), running regressiontests from build system is not supported. Please run gmxtest.pl directly.")
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
    if(DEFINED ENV{CI_JOB_ID} AND GMX_OPENMP)
        set(GMX_REGRESSIONTEST_OMP_THREADS 2 CACHE STRING
            "Number of OpenMP threads to use with the regression tests")
        mark_as_advanced(GMX_REGRESSIONTEST_OMP_THREADS)
        list(APPEND ARGS -ntomp ${GMX_REGRESSIONTEST_OMP_THREADS})
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
    #given that we know they are there and that mdrun might not be executable
    #(e.g. Cray) we enable it.
    list(APPEND ARGS -crosscompile)

    set(REGRESSIONTEST_EXTRA_ARGS "" CACHE STRING 
        "Extra arguments passed to gmxtest")
    mark_as_advanced(REGRESSIONTEST_EXTRA_ARGS)
    list(APPEND ARGS ${REGRESSIONTEST_EXTRA_ARGS})

    if(DEFINED ENV{CI_JOB_ID} AND GMX_THREAD_MPI)
        set(GMX_REGRESSIONTEST_THREAD_MPI_THREADS 2 CACHE STRING
        "Number of thread MPI threads to use for testing")
        mark_as_advanced(GMX_REGRESSIONTEST_THREAD_MPI_THREADS)
        list(APPEND ARGS -nt ${GMX_REGRESSIONTEST_THREAD_MPI_THREADS})
    endif()
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
    foreach(subtest complex essentialdynamics)
        add_test(NAME regressiontests/${subtest}
            #windows requires the command to be perl and not the script
            COMMAND perl "${REGRESSIONTEST_PATH}/gmxtest.pl" ${subtest} ${ARGS})
        set_tests_properties(regressiontests/${subtest} PROPERTIES
            ENVIRONMENT "PATH=${PATH}"
            LABELS "SlowGpuTest"
            RUN_SERIAL TRUE)
    endforeach()
else()
    gmx_add_missing_tests_notice("Regression tests have not been run. If you want to run them from the build system, get the correct version of the regression tests package and set REGRESSIONTEST_PATH in CMake to point to it, or set REGRESSIONTEST_DOWNLOAD=ON.")
endif()


#
# Physical validation tests are opt-in via -DGMX_PHYSICAL_VALIDATION=ON
#
if(GMX_PHYSICAL_VALIDATION)
    include(FindPythonModule)
    # physical validation suite is distributed with the source
    set(PHYSVALTEST_SOURCE_PATH
            "${CMAKE_CURRENT_SOURCE_DIR}/physicalvalidation")
    #       CACHE PATH "Source directory containing physical validation tests.")
    if(NOT EXISTS ${PHYSVALTEST_SOURCE_PATH}/gmx_physicalvalidation.py)
        message(FATAL_ERROR
            "GMX_PHYSICAL_VALIDATION set, but physical validation script not found in ${PHYSVALTEST_SOURCE_PATH}.")
    endif()

    if(CMAKE_CROSSCOMPILING OR CMAKE_CONFIGURATION_TYPES)
        # The following comment is copied from regression tests:
        #     Cross-compiling toolchains require us to compile both front-end and
        #     back-end binaries to run gmxtest.pl.
        # TODO: Look into the details of this.
        # For now, turn it off - our python-gmx interface is probably not that stable for special cases anyway
        message(WARNING
                "With cross-compiling or multi-configuration generators (e.g. Visual Studio),\
                running physical validation tests from build system is not supported.\
                Please run physicalvalidation.py directly.")
        set(GMX_PHYSICAL_VALIDATION OFF CACHE BOOL
                "With cross-compiling or multi-configuration generators, running physical validation tests from build\
                system is not supported." FORCE)

    else()
        #
        # Making sure gmx is built before running physical validation tests
        #
        add_dependencies(run-ctest-phys gmx)
        add_dependencies(run-ctest gmx)

        #
        # Determine arguments passed to physicalvalidation.py
        #
        set(PARGS "")
        list(APPEND PARGS --wd ${CMAKE_CURRENT_BINARY_DIR}/physicalvalidation)
        list(APPEND PARGS --bindir ${CMAKE_BINARY_DIR}/bin)

        if(GMX_LIB_MPI AND NOT MPIEXEC) # find_program failed
            message(WARNING
                "Please set MPIEXEC. Otherwise mpirun is assumed for running tests.")
        endif()
        if(GMX_LIB_MPI)
            # define number of processors in analogy to regression tests?
            if(MPIEXEC)
                list(APPEND ARGS --mpicmd ${MPIEXEC})
            else()
                list(APPEND ARGS --mpicmd mpirun)
            endif()
        endif()
        if(GMX_LIB_MPI)
            message(FATAL_ERROR
                "Physical validation using MPI not supported.")
        endif()

        if(GMX_BINARY_SUFFIX)
            list(APPEND PARGS --suffix ${GMX_BINARY_SUFFIX})
        endif()

        #
        # Give possibility to add args via cache
        #
        set(PHYSVALTEST_EXTRA_ARGS "" CACHE STRING
            "Extra arguments passed to phystest")
        mark_as_advanced(PHYSVALTEST_EXTRA_ARGS)
        list(APPEND PARGS ${PHYSVALTEST_EXTRA_ARGS})

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

        if (NOT Python3_Interpreter_FOUND)
            message(FATAL_ERROR
                    "Python not found. Physical validation requires python. \
                     Install python, set Python3_ROOT_DIR or PYTHON_EXECUTABLE to a valid location, \
                     or set GMX_PHYSICAL_VALIDATION=OFF to disable the physical validation tests.")
        endif()
        #
        foreach(module numpy scipy pymbar) # add further modules if necessary
            find_python_module(${module})
            string(TOUPPER ${module} module_upper)
            if(NOT PYTHONMODULE_${module_upper})
                message(WARNING
                        "Python module ${module} not found. Physical validation relies on ${module}. Make sure\
                        ${module} can be found by PYTHON_EXECUTABLE, or set GMX_PHYSICAL_VALIDATION=OFF to disable\
                        the physical validation tests.")
            endif()
        endforeach()

        #
        # Hook in our own tests
        # Read them from json file to make every system a separate test
        #
        if (GMX_DOUBLE)
            set(PHYSVALTEST_JSON "${PHYSVALTEST_SOURCE_PATH}/systems_d.json")
        else()
            set(PHYSVALTEST_JSON "${PHYSVALTEST_SOURCE_PATH}/systems.json")
        endif()
        file(STRINGS "${PHYSVALTEST_JSON}" json)
        string(REPLACE "\"" "" json ${json})
        string(REPLACE "," "" json ${json})
        string(REPLACE " " "" json ${json})
        string(REPLACE ";" "<<>>" json ${json})
        string(REPLACE "[" "" json ${json})
        string(REPLACE "]" "" json ${json})
        string(REPLACE "{" "" json ${json})
        string(REPLACE "}" "" json ${json})
        string(REPLACE "<<>>" ";" json ${json})
        foreach(line ${json})
            if("${line}" MATCHES "name:")
                string(REPLACE "name:" "" testname ${line})
                add_test(NAME physicalvalidationtests/${testname}
                        COMMAND ${PYTHON_EXECUTABLE} "${PHYSVALTEST_SOURCE_PATH}/gmx_physicalvalidation.py" "${PHYSVALTEST_JSON}" -s ${testname} -a ${PARGS})
                set_tests_properties(physicalvalidationtests/${testname} PROPERTIES
                        ENVIRONMENT "PATH=${PATH}"
                        LABELS "PhysicalValidationTest"
                        RUN_SERIAL TRUE)
            endif()
        endforeach()

        #
        # Create prepare and run targets while all variables are set
        # Will be referenced in CheckTarget.cmake
        #
        # "check-phys-prepare" prepares the systems needed for physical validation for external running
        add_custom_target(check-phys-prepare
                          COMMAND ${PYTHON_EXECUTABLE} "${PHYSVALTEST_SOURCE_PATH}/gmx_physicalvalidation.py" "${PHYSVALTEST_JSON}" -p ${PARGS}
                          COMMENT "Preparing systems for physical validation"
                          DEPENDS gmx)
        # "run-physval-sims" prepares and runs the systems needed for physical validation
        add_custom_target(run-physval-sims
                          COMMAND ${PYTHON_EXECUTABLE} "${PHYSVALTEST_SOURCE_PATH}/gmx_physicalvalidation.py" "${PHYSVALTEST_JSON}" -r ${PARGS}
                          COMMENT "Preparing and running systems for physical validation"
                          DEPENDS gmx)
    endif()
else()
    #
    # Create dummy prepare and run targets
    # Will be referenced in CheckTarget.cmake
    #
    # "check-phys-prepare" prepares the systems needed for physical validation for external running
    add_custom_target(check-phys-prepare
                      COMMAND ${CMAKE_COMMAND} -E echo "NOTE: You called the target `check-phys-prepare`, but ran cmake with\
 `-DGMX_PHYSICAL_VALIDATION=OFF`. The physical validation tests are therefore unavailable,\
 and this target is not doing anything."
                      COMMENT "No physical validation" VERBATIM)
    # "run-physval-sims" prepares and runs the systems needed for physical validation
    add_custom_target(run-physval-sims
                      COMMAND ${CMAKE_COMMAND} -E echo "NOTE: You called the target `run-physval-sims`, but ran cmake with\
 `-DGMX_PHYSICAL_VALIDATION=OFF`. The physical validation tests are therefore unavailable,\
 and this target is not doing anything."
                      COMMENT "No physical validation" VERBATIM)
endif()

gmx_create_missing_tests_notice_target()
