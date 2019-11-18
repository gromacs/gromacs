#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

# Sets version information variables and provides CMake functions for
# generating files based on them
#
# This script provides the following basic version variables that need to be
# maintained manually:
#   GMX_VERSION_MAJOR      Major version number.
#   GMX_VERSION_PATCH      Patch version number.
#       Should always be defined: zero for, e.g., 2016.
#   GMX_VERSION_SUFFIX     String suffix to add to numeric version string.
#       "-dev" is automatically added when not building from a source package,
#       and does not need to be kept here. This mechanism is not quite enough
#       for building a tarball, but setting the CMake cache variable
#       GMX_BUILD_TARBALL=on will suppress the addition of "-dev" to the
#       version string.
#   LIBRARY_SOVERSION_MAJOR so major version for the built libraries.
#       Should be increased for each binary incompatible release. In GROMACS,
#       the typical policy is to increase it at the start of the development
#       cycle for each major/minor version change, but not for patch releases,
#       even if the latter may not always be fully binary compatible.
#       Table of historical values
#         GROMACS     5.0    0
#         GROMACS     5.1    1
#         GROMACS     2016   2
#         GROMACS     2018   3
#         GROMACS     2019   4
#         GROMACS     2020   5
#         GROMACS     2021   6
#   LIBRARY_SOVERSION_MINOR so minor version for the built libraries.
#       Should be increased for each release that changes only the implementation.
#       In GROMACS, the typical policy is to increase it for each patch version
#       change, even if they may not always be fully binary compatible.
#       If it is somehow clear that the ABI implementation has not changed
#       in a patch release, this variable should not increase. Release candidate
#       and beta versions will not increase this number, since nobody should
#       write code against such versions.
#   LIBRARY_VERSION        Full library version.
#   REGRESSIONTEST_BRANCH  For builds not from source packages, name of the
#       regressiontests branch at gerrit.gromacs.org whose HEAD can test this
#       code, *if* this code is recent enough (i.e., contains all changes from
#       the corresponding code branch that affects the regression test
#       results). Even after a release branch is forked for the source
#       repository, the correct regressiontests branch can still be master,
#       because we do not fork it until behaviour needs to change.
#   REGRESSIONTEST_MD5SUM
#       The MD5 checksum of the regressiontest tarball. Only used when building
#       from a source package.
#   GMX_SOURCE_DOI_ID
#       ID collected from Zenodo connected to the doi for a released version
#       used to identify the source when building an official released version.
#       This ID is used for the source code tarball.
#   GMX_MANUAL_DOI_ID
#       Same as above, but for the reference manual.
# They are collected into a single section below.
# The following variables are set based on these:
#   GMX_VERSION            String composed from GMX_VERSION_* numeric variables
#       above. Example: 4.6.1, 5.0, 2016
#   GMX_VERSION_STRING     String with GMX_VERSION suffixed with the given
#       suffix and possibly "-dev" for builds not from a source package.
#   GMX_VERSION_NUMERIC    Numeric version number (e.g., 40601 for 4.6.1, 20160001 for 2016.1).
#   GMX_API_VERSION        Numeric API version.
#       This is currently set automatically to GMX_VERSION_NUMERIC, but may
#       become manually maintained in the future if there will be releases
#       where the API does not change, but programs/libraries do.
#       In such a case, this should be the first version where the current API
#       appeared.
#   REGRESSIONTEST_VERSION For source packages, version number of the
#       matching regressiontests tarball.  Not used for builds not from source
#       packages.
# The latter two are used to generate gromacs/version.h to allow software
# written against the GROMACS API to provide some #ifdef'ed code to support
# multiple GROMACS versions.
#
# This script also declares machinery to generate and obtain version
# information from a git repository.  This is enabled by default if the source
# tree is a git, but can be disabled with
#   GMX_GIT_VERSION_INFO           Advanced CMake variable to disable git
#                                  version info generation.
# If the version generation is disabled, then the source and manual doi
# will be based on the stored values for the ID.
# The main interface to this machinery is the gmx_configure_version_file()
# CMake function.  The signature is
#   gmx_configure_version_file(<input> <output>
#                              [REMOTE_HASH]
#                              [TARGET <target>]
#                              [COMMENT <comment>])
#   <input>      Specify the input and output files as for configure_file().
#   <output>     The configuration is done with configure_file(... @ONLY) with
#                the following variables defined (as well as all the
#                GMX_VERSION* variables from above):
#                  GMX_VERSION_STRING_FULL
#                  GMX_VERSION_FULL_HASH
#                  GMX_VERSION_CENTRAL_BASE_HASH
#                The output file is created during build time, so any dependent
#                targets should specify it as a dependency.
#   REMOTE_HASH  Currently, this has no effect, but it signifies that the
#                <input> file is using the CENTRAL_BASE_HASH variable.
#                This variable is much more expensive to initialize than the
#                others, so this allows local changes in this file to only
#                compute that value when required if that becomes necessary.
#   TARGET       By default, this function uses add_custom_command() to
#                generate the output file.  If TARGET is specified, then
#                add_custom_target() is used to create a target with the given
#                name <target> that runs this custom command.  Use this if
#                the same file will be used for multiple downstream targets,
#                or if the explicit target for the file is otherwise
#                necessary.
#   COMMENT      Set a custom comment to be shown when building the rule
#                (see add_custom_command(... COMMENT <comment>)).
# As an alternative to using this script, also the following variables are
# provided (can be useful when generating more complex CMake scripts that do
# build-time tasks):
#   VERSION_INFO_CMAKE_SCRIPT
#       Absolute path to a CMake script that can be included using include()
#       to declare the GMX_VERSION_* variables documented for
#       gmx_configure_version_file().
#   VERSION_INFO_DEPS
#       If a custom command depends on VERSION_INFO_CMAKE_SCRIPT, then it
#       should add ${VERSION_INFO_DEPS} to its DEPENDS list to get the
#       appropriate dependencies.
# TODO: If someone wants to add a custom target that depends on
# VERSION_INFO_CMAKE_SCRIPT, a separate variable may be needed for those
# dependencies.
#
# The version string printed by 'gmx -version' (and also printed in the startup
# header) can provide useful information for, e.g., diagnosing bug reports and
# identifying what exact version the user was using.  The following formats are
# possible (with examples given for a particular version):
#   2018.1       Plain version number without any suffix signifies a build from
#                a released source tarball.
#   2018.1-dev   '-dev' suffix signifies all other builds. If there is no other
#                information, either the user built the code outside any git
#                repository, or disabled the version info generation.
#   2018.1-dev-YYYYMMDD-1234abc
#                The YYYYMMDD part shows the commit date (not author date) of
#                the HEAD commit from which the code was built.  The abbreviated
#                hash is the hash of that commit (the full hash is available in
#                'gmx -version' output).
#                If the HEAD hash is not identified as coming from branches in
#                "authoritative" GROMACS repositories, 'gmx -version' will show
#                the nearest ancestor commit that is identified as such (but see
#                the '-local' and '-unknown' suffixes below).
#   2018.1-dev-YYYYMMDD-1234abc-dirty
#                As above, but there were local modifications in the source tree
#                when the code was built.
#   2018.1-dev-YYYYMMDD-1234abc-unknown
#                As above, but there were no remotes in the repository that
#                could be identified as "authoritative" GROMACS repositories.
#                This happens if the code is not cloned from git.gromacs.org
#                or gerrit.gromacs.org.
#   2018.1-dev-YYYYMMDD-1234abc-local
#                As above, but there were no commits in the recent history of
#                the branch that could be identified as coming from
#                "authoritative" GROMACS repositories.  This should be
#                relatively rare.
#
# Other variables set here are not intended for use outside this file.
# The scripts gmxGenerateVersionInfo.cmake and gmxConfigureVersionInfo.cmake
# are used internally by this machinery, as well as VersionInfo.cmake.cmakein.

#####################################################################
# Manually maintained version info

# The GROMACS convention is that these are the version number of the next
# release that is going to be made from this branch.
set(GMX_VERSION_MAJOR 2021)
set(GMX_VERSION_PATCH 0)
# The suffix, on the other hand, is used mainly for betas and release
# candidates, where it signifies the most recent such release from
# this branch; it will be empty before the first such release, as well
# as after the final release is out.
set(GMX_VERSION_SUFFIX "")

# Conventionally with libtool, any ABI change must change the major
# version number, the minor version number should change if it's just
# the implementation that has been altered, and the third number
# counts the number of old major versions that will still run if
# linked to this library (i.e. it is not a patch number). See the
# above descriptions of LIBRARY_SOVERSION_* for policy for changes
# here. The important thing is to minimize the chance of third-party
# code being able to dynamically link with a version of libgromacs
# that might not work.
set(LIBRARY_SOVERSION_MAJOR 6)
set(LIBRARY_SOVERSION_MINOR 0)
set(LIBRARY_VERSION ${LIBRARY_SOVERSION_MAJOR}.${LIBRARY_SOVERSION_MINOR}.0)

#####################################################################
# General version management based on manually set numbers

if (GMX_VERSION_PATCH)
    set(GMX_VERSION "${GMX_VERSION_MAJOR}.${GMX_VERSION_PATCH}")
else()
    set(GMX_VERSION "${GMX_VERSION_MAJOR}")
endif()
set(GMX_VERSION_STRING "${GMX_VERSION}${GMX_VERSION_SUFFIX}")

# If you are making a custom fork of GROMACS, please describe your
# fork, perhaps with its version number, in the value of
# GMX_VERSION_STRING_OF_FORK here. This string will appear in the
# header of log files that mdrun writes. This will help you, your
# users, your system administrators, your maintainers and the
# maintainers of GROMACS core understand how to troubleshoot and
# reproduce potential problems.
#
# If you are distributing a patch to GROMACS, then this change would
# be great as part of your patch. Otherwise for personal use, you can
# also just set a CMake cache variable.
set(GMX_VERSION_STRING_OF_FORK "" CACHE INTERNAL
    "Version string for forks of GROMACS to set to describe themselves")
mark_as_advanced(GMX_VERSION_STRING_OF_FORK)
if (GMX_VERSION_STRING_OF_FORK)
    set(GMX_VERSION_STRING "${GMX_VERSION_STRING}-${GMX_VERSION_STRING_OF_FORK}")
endif()

option(GMX_BUILD_TARBALL "Build tarball without -dev version suffix" OFF)
mark_as_advanced(GMX_BUILD_TARBALL)
# If run with cmake -P, the -dev suffix is managed elsewhere.
if (NOT SOURCE_IS_SOURCE_DISTRIBUTION AND
    NOT GMX_BUILD_TARBALL AND
    NOT CMAKE_SCRIPT_MODE_FILE)
    set(GMX_VERSION_STRING "${GMX_VERSION_STRING}-dev")
endif()

set(REGRESSIONTEST_VERSION "${GMX_VERSION_STRING}")
set(REGRESSIONTEST_BRANCH "refs/heads/master")
# Run the regressiontests packaging job with the correct pakage
# version string, and the release box checked, in order to have it
# build the regressiontests tarball with all the right naming. The
# naming affects the md5sum that has to go here, and if it isn't right
# release workflow will report a failure.
set(REGRESSIONTEST_MD5SUM "162edf7d9d520672f5fa5888aae60989" CACHE INTERNAL "MD5 sum of the regressiontests tarball for this GROMACS version")

math(EXPR GMX_VERSION_NUMERIC
     "${GMX_VERSION_MAJOR}*10000 + ${GMX_VERSION_PATCH}")
set(GMX_API_VERSION ${GMX_VERSION_NUMERIC})

# If run with cmake -P from releng scripts, print out necessary version info
# as JSON.
if (CMAKE_SCRIPT_MODE_FILE)
    message("{ \"version\": \"${GMX_VERSION_STRING}\", \"regressiontest-md5sum\": \"${REGRESSIONTEST_MD5SUM}\" }")
    return()
endif()

# Set those values only in release versions, after getting the identifiers
# from Zenodo for the manual and source code
# Has to be done by hand before every final release
# Use force to override anything given as a cmake command line input
# Actual input depends on the GMX_VERSION_STRING_OF_FORK variable being set or not.
# If it is set, we always default to an empty string, otherwise to the value set for the release build.
if (GMX_VERSION_STRING_OF_FORK)
    set(GMX_MANUAL_DOI_INTERNAL "")
    set(GMX_SOURCE_DOI_INTERNAL "")
else()
    set(GMX_MANUAL_DOI_INTERNAL "") # Set correct doi string here
    set(GMX_SOURCE_DOI_INTERNAL "") # Set correct doi string here
endif()
set(GMX_MANUAL_DOI ${GMX_MANUAL_DOI_INTERNAL} CACHE INTERNAL "reserved doi for GROMACS manual" FORCE)
set(GMX_SOURCE_DOI ${GMX_SOURCE_DOI_INTERNAL} CACHE INTERNAL "reserved doi for GROMACS source code" FORCE)

#####################################################################
# git version info management

# There can be clusters where git and CMake can run on nodes where the other is
# not available, accessing the same source tree.
# Should be unlikely, but doesn't hurt to check.
set(_git_info_default OFF)
if (SOURCE_IS_GIT_REPOSITORY)
    find_package(Git)
    if (GIT_FOUND)
        set(_git_info_default ON)
    endif()
endif()
option(GMX_GIT_VERSION_INFO "Generate git version information" ${_git_info_default})
mark_as_advanced(GMX_GIT_VERSION_INFO)
# Detect preconditions for version info generation if it is requested.
if (GMX_GIT_VERSION_INFO)
    if (NOT SOURCE_IS_GIT_REPOSITORY)
        message(FATAL_ERROR
            "Cannot generate git version information from source tree not under git. "
            "Set GMX_GIT_VERSION_INFO=OFF to proceed.")
    endif()
    # We need at least git v1.5.3 be able to parse git's date output.
    if (NOT GIT_FOUND OR GIT_VERSION_STRING VERSION_LESS "1.5.3")
        message(FATAL_ERROR
            "No compatible git version found (>= 1.5.3 required). "
            "Won't be able to generate development version information. "
            "Set GMX_GIT_VERSION_INFO=OFF to proceed.")
    endif()
endif()

include(gmxCustomCommandUtilities)

# The first two are also for use outside this file, encapsulating the details
# of how to use the generated VersionInfo.cmake.
set(VERSION_INFO_CMAKE_FILE   ${PROJECT_BINARY_DIR}/VersionInfo.cmake)
set(VERSION_INFO_DEPS         ${VERSION_INFO_CMAKE_FILE})
# Capture the location of the necessary files in internal variables for use in
# the function below.
set(VERSION_INFO_CMAKEIN_FILE     ${CMAKE_CURRENT_LIST_DIR}/VersionInfo.cmake.cmakein)
set(VERSION_INFO_CONFIGURE_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/gmxConfigureVersionInfo.cmake)
# A set of directories to scan for calculating the hash of source files.
set(SET_OF_DIRECTORIES_TO_CHECKSUM  "${PROJECT_SOURCE_DIR}/src")
list(APPEND SET_OF_DIRECTORIES_TO_CHECKSUM "${PROJECT_SOURCE_DIR}/python_packaging")
# Try to find python for the checksumming script
set(PythonInterp_FIND_QUIETLY ON)
find_package(PythonInterp 3.5)

# Rules to create the VersionInfo.cmake file.
# For git info, the sequence is:
#   1. (configure time) VersionInfo.cmake.cmakein -> VersionInfo-partial.cmake.cmakein
#        - Set all variables that are known at configure time.
#   2. (build time)     VersionInfo-partial.cmake.cmakein -> VersionInfo.cmake
#        - Set variables that may change as result of repository state changes
#          (i.e., everything that requires running git).
#        - Runs every time as a git-version-info target, but the output file
#          timestamp only changes if its contents actually change.
#        - Depending on the native build system, this may run once per build
#          or once per each time it is required for step 3.
#   3. (build time)     VersionInfo.cmake -> other files
#        - Set variables in files specified with gmx_configure_version_file()
#          using the values generated in step 2.
#        - Each file runs as a custom command that depends on the previous
#          steps, and runs only if the VersionInfo.cmake file is newer than the
#          output file.
# Without git info, the sequence is:
#  1. (configure time) VersionInfo.cmake.cmakein -> VersionInfo.cmake
#        - Everything is known at configure time, so the output is generated
#          immediately with all variables set (git info will be empty).
#  2. (build time)     VersionInfo.cmake -> other files
#        - As with git info, processes files from gmx_configure_version_file().
#        - These are again custom commands that depend on the output from
#          step 1, so they get regenerated only when the static version info
#          changes.
if (GMX_GIT_VERSION_INFO)
    # Configure information known at this time into a partially filled
    # version info file.
    set(VERSION_INFO_CMAKEIN_FILE_PARTIAL
        ${PROJECT_BINARY_DIR}/VersionInfo-partial.cmake.cmakein)
    # Leave these to be substituted by the custom target below.
    set(GMX_VERSION_STRING_FULL       "\@GMX_VERSION_STRING_FULL\@")
    set(GMX_VERSION_FULL_HASH         "\@GMX_VERSION_FULL_HASH\@")
    set(GMX_VERSION_CENTRAL_BASE_HASH "\@GMX_VERSION_CENTRAL_BASE_HASH\@")
    configure_file(${VERSION_INFO_CMAKEIN_FILE}
                   ${VERSION_INFO_CMAKEIN_FILE_PARTIAL}
                   @ONLY)
    # If generating the version info, create a target that runs on every build
    # and does the actual git calls, storing the results into a CMake script.
    # This needs to be run at build time to update the version information
    # properly when the git hash changes, but the build system does not.
    # All targets added by gmx_configure_version_file() use the information
    # from this script to get their variables from, removing the need to run
    # git multiple times and simplifying reuse for other purposes.
    gmx_add_custom_output_target(git-version-info RUN_ALWAYS
        OUTPUT ${VERSION_INFO_CMAKE_FILE}
        COMMAND ${CMAKE_COMMAND}
            -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
            -D PROJECT_VERSION=${GMX_VERSION_STRING}
            -D PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
            -D VERSION_CMAKEIN=${VERSION_INFO_CMAKEIN_FILE_PARTIAL}
            -D VERSION_OUT=${VERSION_INFO_CMAKE_FILE}
            -P ${CMAKE_CURRENT_LIST_DIR}/gmxGenerateVersionInfo.cmake
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating git version information")
    list(APPEND VERSION_INFO_DEPS git-version-info)
else()
    # If the version info is static, just generate the CMake script with the
    # version variables during the CMake run.
    set(GMX_VERSION_STRING_FULL       ${GMX_VERSION_STRING})
    set(GMX_VERSION_FULL_HASH         "")
    set(GMX_VERSION_CENTRAL_BASE_HASH "")
    # To notify the user during compilation and at runtime that the build source
    # has not been modified after unpacking the source tarball, the contents are hashed
    # to be compared to a hash computed during the release process. If the hash matches
    # all is fine and the user gets a message in the log file indicating that.
    # If either the release hash file is missing, or if the hash does not match
    # a different message is printed to indicate that the source has been changed
    # compared to the version actually released. This is not needed in case a build
    # is done in git, as we have the information there already.
    # This is not done if the user has explicitly set an additional custom version string with
    # -DGMX_VERSION_STRING_OF_FORK, as this indicates that they are knowing that a custom
    # version of GROMACS is in use.
    set(RELEASE_CHECKSUM_FILE "${PROJECT_SOURCE_DIR}/src/reference_checksum")
    if(NOT GMX_VERSION_STRING_OF_FORK OR "${GMX_VERSION_STRING_OF_FORK}" STREQUAL "")
        if(EXISTS ${RELEASE_CHECKSUM_FILE} AND PythonInterp_FOUND)
            file(READ ${RELEASE_CHECKSUM_FILE} GMX_RELEASE_SOURCE_FILE_CHECKSUM)
            string(STRIP ${GMX_RELEASE_SOURCE_FILE_CHECKSUM} GMX_RELEASE_SOURCE_FILE_CHECKSUM)
            set(CHECKSUM_RESULT_FILE "${CMAKE_CURRENT_BINARY_DIR}/computed_checksum")
            execute_process(COMMAND ${PYTHON_EXECUTABLE} 
                                    ${PROJECT_SOURCE_DIR}/admin/createFileHash.py
                                    -s ${SET_OF_DIRECTORIES_TO_CHECKSUM}
                                    -o ${CHECKSUM_RESULT_FILE}
                            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                            OUTPUT_QUIET)
                        file(READ ${CHECKSUM_RESULT_FILE} GMX_CURRENT_SOURCE_FILE_CHECKSUM)
            string(STRIP ${GMX_CURRENT_SOURCE_FILE_CHECKSUM} GMX_CURRENT_SOURCE_FILE_CHECKSUM)
            if(NOT ${GMX_RELEASE_SOURCE_FILE_CHECKSUM} STREQUAL ${GMX_CURRENT_SOURCE_FILE_CHECKSUM})
                set(GMX_VERSION_STRING_FULL "${GMX_VERSION_STRING_FULL}-MODIFIED")
                message(STATUS "The source code for this GROMACS installation is different from the officially released version.")
            endif()
        elseif(PythonInterp_FOUND)
            set(GMX_VERSION_STRING_FULL "${GMX_VERSION_STRING_FULL}-UNCHECKED")
            set(GMX_RELEASE_SOURCE_FILE_CHECKSUM "NoChecksumFile")
            set(GMX_CURRENT_SOURCE_FILE_CHECKSUM "NoChecksumFile")
            message(WARNING "Could not valdiate the GROMACS source due to missing reference checksum file.")
        else()
            set(GMX_VERSION_STRING_FULL "${GMX_VERSION_STRING_FULL}-UNCHECKED")
            set(GMX_RELEASE_SOURCE_FILE_CHECKSUM "NoPythonAvailable")
            set(GMX_CURRENT_SOURCE_FILE_CHECKSUM "NoPythonAvailable")
            message(STATUS "Could not calculate checksum of source files without Python")
        endif()
    endif()
    configure_file(${VERSION_INFO_CMAKEIN_FILE} ${VERSION_INFO_CMAKE_FILE})
endif()
unset(GMX_VERSION_STRING_FULL)
unset(GMX_VERSION_FULL_HASH)
unset(GMX_VERSION_CENTRAL_BASE_HASH)
unset(GMX_RELEASE_SOURCE_FILE_CHECKSUM)
unset(GMX_CURRENT_SOURCE_FILE_CHECKSUM)


# What file the checksum should be written to
set(CHECKSUM_FILE "${PROJECT_SOURCE_DIR}/src/reference_checksum")

# Target that allows checksumming a source tree when producing a tarball.
# Allows verification of builds from the tarball to make sure the source had
# not been tampered with.
# Note: The RUN_ALWAYS here is to regenerate the hash file only, it does not
# mean that the target is run in all builds
if (PYTHONINTERP_FOUND)
    gmx_add_custom_output_target(reference_checksum RUN_ALWAYS
        OUTPUT ${CHECKSUM_FILE}
        COMMAND ${PYTHON_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/admin/createFileHash.py
            -s ${SET_OF_DIRECTORIES_TO_CHECKSUM}
            -o ${CHECKSUM_FILE}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating reference checksum of source files")
else()
    add_custom_target(reference_checksum
        COMMAND ${CMAKE_COMMAND} -E echo
        "Can not checksum files without python being available"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating reference checksum of source files")
endif()

# The main user-visible interface to the machinery.
# See documentation at the top of the script.
function (gmx_configure_version_file INFILE OUTFILE)
    include(CMakeParseArguments)
    set(_options REMOTE_HASH)
    set(_one_value_args COMMENT TARGET)
    set(_multi_value_args EXTRA_VARS)
    cmake_parse_arguments(
        ARG "${_options}" "${_one_value_args}" "${_multi_value_args}" ${ARGN})
    if (ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unknown arguments: ${ARG_UNPARSED_ARGUMENTS}")
    endif()
    # Some callers may pass partial paths that do not really make sense,
    # so create a default comment that only contains the actual file name.
    get_filename_component(_basename ${OUTFILE} NAME)
    set(_comment "Generating ${_basename}")
    if (ARG_COMMENT)
        set(_comment ${ARG_COMMENT})
    endif()
    # Mimic configure_file()
    if (NOT IS_ABSOLUTE ${INFILE})
        set(INFILE ${CMAKE_CURRENT_SOURCE_DIR}/${INFILE})
    endif()
    # Create command-line definitions for the requested variables
    set(_extra_var_defines)
    foreach(_var ${ARG_EXTRA_VARS})
        list(APPEND _extra_var_defines -D "${_var}=${${_var}}")
    endforeach()
    # The touch command is necessary to ensure that after the target is run,
    # the timestamp is newer than in the input files.
    add_custom_command(OUTPUT ${OUTFILE}
        COMMAND ${CMAKE_COMMAND}
            -D VERSION_VARIABLES=${VERSION_INFO_CMAKE_FILE}
            -D VERSION_CMAKEIN=${INFILE}
            -D VERSION_OUT=${OUTFILE}
            ${_extra_var_defines}
            -P ${VERSION_INFO_CONFIGURE_SCRIPT}
        COMMAND ${CMAKE_COMMAND} -E touch ${OUTFILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${INFILE} ${VERSION_INFO_DEPS} ${VERSION_INFO_CONFIGURE_SCRIPT}
        COMMENT "${_comment}"
        VERBATIM)
    if (ARG_TARGET)
        add_custom_target(${ARG_TARGET} DEPENDS ${OUTFILE} VERBATIM)
        gmx_set_custom_target_output(${ARG_TARGET} ${OUTFILE})
    endif()
endfunction()
