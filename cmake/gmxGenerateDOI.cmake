#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

# Generate Gromacs DOI for new releases.

# This script will automatically register the new software version
# tarball and manual with a Zenodo doi upon release, so that it can
# be properly referenced in the future.

# This is mainly done by calling the python scripts
# reserve_doi.py and publish_doi.py at different stages of the
# build. The doi string is included in the tarball sources
# and the manual for this release.

#
# Paul Bauer (paul.bauer.q@gmail.com)

# function that actually gets the values from Zenodo
function (gmx_register_doi REGTYPE)
    include(CMakeParseArguments)
    set(_options REMOTE_HASH)
    set(_one_value_args COMMENT TARGET)
    set(_multi_value_args EXTRA_VARS)
    cmake_parse_arguments(
        ARG "${_options}" "${_one_value_args}" "${_multi_value_args}" ${ARGN})
    if (ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unknown arguments: ${ARG_UNPARSED_ARGUMENTS}")
    endif()
    # get doi strings for source tarball and manual
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/admin/reserve_doi.py ${REGTYPE} ${GMX_VERSION_STRING} ${TOKEN_PATH}
        OUTPUT_VARIABLE GMX_DOI_OUTPUT_VARIABLE
        RESULT_VARIABLE GMX_DOI_RESULT_VARIABLE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    if(NOT ${GMX_DOI_RESULT_VARIABLE} EQUAL "0")
    # return with error here
      MESSAGE(FATAL_ERROR "Something went wrong while registering the doi for the ${REGTYPE}.\nOutput value is: ${GMX_DOI_OUTPUT_VARIABLE}")
    endif()
    # need to split return string into two parts
    separate_arguments(GMX_DOI_OUTPUT_VARIABLE)
    set(RETURNVAR)
    string(TOUPPER ${REGTYPE} REGTYPE)
    # first element of list is always the doi string
    list(GET GMX_DOI_OUTPUT_VARIABLE 0 RETURNVAR)
    set(GMX_${REGTYPE}_DOI ${RETURNVAR} CACHE INTERNAL "reserved doi for GROMACS ${REGTYPE}" )
    # second list element is the submission id
    list(GET GMX_DOI_OUTPUT_VARIABLE 1 RETURNVAR)
    set(GMX_${REGTYPE}_ID ${RETURNVAR} CACHE INTERNAL "reserved doi id for GROMACS ${REGTYPE}")
endfunction()

# function that allows the interface to the update_doi.py python script
# to update entries on zenodo
function (gmx_update_doi REGTYPE OLD_ID)
    include(CMakeParseArguments)
    set(_options REMOTE_HASH)
    set(_one_value_args COMMENT TARGET)
    set(_multi_value_args EXTRA_VARS)
    cmake_parse_arguments(
        ARG "${_options}" "${_one_value_args}" "${_multi_value_args}" ${ARGN})
    if (ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unknown arguments: ${ARG_UNPARSED_ARGUMENTS}")
    endif()
    # get doi strings for source tarball and manual
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/admin/update_doi.py ${REGTYPE} ${GMX_VERSION_STRING} ${OLD_ID} ${TOKEN_PATH}
        OUTPUT_VARIABLE GMX_DOI_OUTPUT_VARIABLE
        RESULT_VARIABLE GMX_DOI_RESULT_VARIABLE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    if(NOT ${GMX_DOI_RESULT_VARIABLE} EQUAL "0")
    # return with error here
      MESSAGE(FATAL_ERROR "Something went wrong while registering the doi for the ${REGTYPE}.\nOutput value is: ${GMX_DOI_OUTPUT_VARIABLE}")
    endif()
    # need to split return string into two parts
    separate_arguments(GMX_DOI_OUTPUT_VARIABLE)
    set(RETURNVAR)
    string(TOUPPER ${REGTYPE} REGTYPE)
    # first element of list is always the doi string
    list(GET GMX_DOI_OUTPUT_VARIABLE 0 RETURNVAR)
    set(GMX_${REGTYPE}_DOI ${RETURNVAR} CACHE INTERNAL "reserved doi for GROMACS ${REGTYPE}" )
    # second list element is the submission id
    list(GET GMX_DOI_OUTPUT_VARIABLE 1 RETURNVAR)
    set(GMX_${REGTYPE}_ID ${RETURNVAR} CACHE INTERNAL "reserved doi id for GROMACS ${REGTYPE}")
endfunction()

# master variables controlling build and submission
option(GMX_GET_RELEASE_DOI "Activate for a build of the release tarball to submit for publishing using Zenodo" OFF)
mark_as_advanced(GMX_GET_RELEASE_DOI)

# set file path variables so they are available everywhere below
set(INFILE "${CMAKE_SOURCE_DIR}/cmake/gmxDOIVersion.cmake.cmakein")
set(OUTFILE "${CMAKE_SOURCE_DIR}/cmake/gmxDOIVersion.cmake")

# Only do this if we are building the release binary
# Instead, check for existing gmxDOIVersion.cmake file and use it to populate
# the variables for source and manual, if build from the source code tarball
# if the file does not exist, the build is either from a git repository or
# an archive not generated from a release build. In this case, the file
# is generated with empty variables
if(NOT GMX_GET_RELEASE_DOI)
    # check if the populated gmxDOIVersion file exists or not
    # if not, we are not originating from a release build and the 
    # doi strings should be empty
    if(NOT EXISTS ${OUTFILE})
        set(GMX_MANUAL_DOI "" CACHE INTERNAL "reserved doi for GROMACS manual")
        set(GMX_SOURCE_DOI "${GMX_MANUAL_DOI}" CACHE INTERNAL "reserved doi for GROMACS manual")
    else()
        # Source code archive that originated from a release build.
        # gmxDOIVersion.cmake is generated from release builds,
        # and populates both manual and source doi strings. 
        include(gmxDOIVersion)
        endif()
        add_custom_target(gmx-publish-source
            COMMENT "Dummy target"
            )
        add_custom_target(gmx-publish-manual
            COMMENT "Dummy target"
            )
else()
    # bool values to control if manual needs to be registered again or not
    set(GMX_HAS_MANUAL_DOI FALSE)
    set(GMX_HAS_SOURCE_DOI FALSE)
    if(DEFINED GMX_MANUAL_DOI)
        set(GMX_HAS_MANUAL_DOI TRUE)
        # Manual already registered, skipping
    endif()
    if(DEFINED GMX_SOURCE_DOI)
        set(GMX_HAS_SOURCE_DOI TRUE)
        # Source code archive already registered, skipping
    endif()
    # make sure we have python available
    find_package(PythonInterp REQUIRED)
    # The Source_package_master build is configured to make available
    # to the build job a secret file (from the Jenkins credential file
    # gmxtoken.txt), and to put the name of the actual file in the
    # environment variable TOKEN_PATH.
    set(TOKEN_PATH $ENV{ZenodoTokenFile})


    if(GMX_UPDATE_RELEASE)
        # we want to update a previous release, so we need to run a different script
        # this will populate the remaining variables as before, but will use the old information
        # to set stuff such as metadata
        # only supported for source code updates, manual should stay the same, and gets its 
        # doi information from the same script through a shorter code path for inclusion in the files
        # but no new target is registered to publish the new version

        # set id values for previous versions of manual and source code from
        # environment variables InputManualID and InputSourceID
        if(NOT ${GMX_HAS_MANUAL_DOI})
            set(GMX_REGTYPE "manual")
            set(INPUT_ID  $ENV{InputManualID})
            gmx_update_doi(${GMX_REGTYPE} ${INPUT_ID})
        endif()
        if(NOT ${GMX_HAS_SOURCE_DOI})
            set(GMX_REGTYPE "source")
            set(INPUT_ID  $ENV{InputSourceID})
            gmx_update_doi(${GMX_REGTYPE} ${INPUT_ID})
        endif()
    else()
        # only run if we did not populate the variables in a previous build
        # unlikely to happen, but better be safe than sorry
        if(NOT ${GMX_HAS_MANUAL_DOI})
            set(GMX_REGTYPE "manual")
            gmx_register_doi(${GMX_REGTYPE})
        endif()
        # see above, not overwriting exisiting variables
        if(NOT ${GMX_HAS_SOURCE_DOI})
            set(GMX_REGTYPE "source")
            gmx_register_doi(${GMX_REGTYPE})
        endif()
    endif()
    # populate gmxDOIVersion file so that doi string is persistent
    # in the release source code archive
    configure_file(${INFILE} ${OUTFILE})
    add_custom_target(gmx-publish-source
        COMMAND
            ${PYTHON_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/admin/publish_doi.py
            source
            ${GMX_VERSION_STRING}
            ${GMX_SOURCE_ID}
            ${CMAKE_CURRENT_BINARY_DIR}
            ${TOKEN_PATH}
        COMMENT "Final publishing of GROMACS source done with Zenodo"
        VERBATIM
        )
    if (NOT GMX_UPDATE_RELEASE)
    add_custom_target(gmx-publish-manual
        COMMAND
            ${PYTHON_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/admin/publish_doi.py
            manual
            ${GMX_VERSION_STRING}
            ${GMX_MANUAL_ID}
            ${CMAKE_CURRENT_BINARY_DIR}
            ${TOKEN_PATH}
        COMMENT "Final publishing of GROMACS manual done with Zenodo"
        VERBATIM
        )
    else()
        add_custom_target(gmx-publish-manual
            COMMENT "Dummy target"              
            )                                   
    endif()
endif()
