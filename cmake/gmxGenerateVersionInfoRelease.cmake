#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020, by the GROMACS development team, led by
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

# Generate GROMACS release build version information, as well as build
# directory integrity checking.

# This script generates version information for a build from a release tarball
# source tree based on the release information in the released tarball.
# It is assumed that by default the script is run in cmake script mode.
# If *not* called in script mode but used in generating cache variables,
# GEN_VERSION_INFO_INTERNAL has to be set ON.
#
# The following variables have to be previously defined:
# PROJECT_VERSION         - hard-coded version string (generated info is appended)
# PROJECT_SOURCE_DIR      - top level source directory
# DIRECTORIES_TO_CHECKSUM - List of directories to hash
# VERSION_CMAKEIN         - path to an input template file
# VERSION_OUT             - path to the output file
#
# The following a possible additional definitions
# PYTHON_EXECUTABLE       - Needed to run checking stage of current tree
# VERSION_STRING_OF_FORK  - Possibly defined custom version string to override
#                          process of checking source file hashes.
# Output:
# VERSION_OUT is configured from the input VERSION_CMAKEIN
# using the variables listed below.
#
# GMX_VERSION_STRING_FULL       - version string
# GMX_RELEASE_SOURCE_FILE_CHECKSUM - Sha256 hash of source tree files at release
# GMX_CURRENT_SOURCE_FILE_CHECKSUM - Sha256 hash of current source tree files
#
# Paul Bauer (paul.bauer.q@gmail.com)
# (authors of git Version of the script that this is based on below)
# Szilard Pall (pszilard@cbr.su.se)
# Teemu Murtola (teemu.murtola@gmail.com)

# Check input variables.
if("${PROJECT_VERSION}" STREQUAL "")
    message(FATAL_ERROR "PROJECT_VERSION undefined!")
endif()
if (NOT EXISTS "${PROJECT_SOURCE_DIR}")
    message(FATAL_ERROR "Project source directory ${PROJECT_SOURCE_DIR} does not exist")
endif()
if ("${DIRECTORIES_TO_CHECKSUM}" STREQUAL "")
    message(FATAL_ERROR "Need list of directories to generate the hash for")
endif()
if ("${VERSION_CMAKEIN}" STREQUAL "")
    message(FATAL_ERROR "Missing input parameter VERSION_CMAKEIN!")
endif()
if ("${VERSION_OUT}" STREQUAL "")
    message(FATAL_ERROR "Missing input parameter VERSION_OUT!")
endif()

# Prepare version string to populate
set(GMX_VERSION_STRING_FULL          ${PROJECT_VERSION})

# We had to pass the directory list as a string, so now we convert it back to a list
string(REPLACE ":" ";" DIRECTORIES_TO_CHECKSUM_LIST ${DIRECTORIES_TO_CHECKSUM})

# Prepare for checking source tree file hashes.
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
if(NOT VERSION_STRING_OF_FORK OR "${VERSION_STRING_OF_FORK}" STREQUAL "")
    if(EXISTS ${RELEASE_CHECKSUM_FILE} AND PYTHON_EXECUTABLE)
        file(READ ${RELEASE_CHECKSUM_FILE} GMX_RELEASE_SOURCE_FILE_CHECKSUM)
        string(STRIP ${GMX_RELEASE_SOURCE_FILE_CHECKSUM} GMX_RELEASE_SOURCE_FILE_CHECKSUM)
        set(CHECKSUM_RESULT_FILE "${CMAKE_CURRENT_BINARY_DIR}/computed_checksum")
        execute_process(COMMAND ${PYTHON_EXECUTABLE}
                                ${PROJECT_SOURCE_DIR}/admin/createFileHash.py
                                -s ${DIRECTORIES_TO_CHECKSUM_LIST}
                                -o ${CHECKSUM_RESULT_FILE}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                        OUTPUT_QUIET)
                    file(READ ${CHECKSUM_RESULT_FILE} GMX_CURRENT_SOURCE_FILE_CHECKSUM)
        string(STRIP ${GMX_CURRENT_SOURCE_FILE_CHECKSUM} GMX_CURRENT_SOURCE_FILE_CHECKSUM)
        if(NOT ${GMX_RELEASE_SOURCE_FILE_CHECKSUM} STREQUAL ${GMX_CURRENT_SOURCE_FILE_CHECKSUM})
            set(GMX_VERSION_STRING_FULL "${GMX_VERSION_STRING_FULL}-MODIFIED")
            message(STATUS "The source code for this GROMACS installation is different from the officially released version.")
        endif()
    elseif(PYTHON_EXECUTABLE)
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

# Generate the output file.
configure_file(${VERSION_CMAKEIN} ${VERSION_OUT})
