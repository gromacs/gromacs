#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

# Generate Gromacs development build version information.

# This script generates version information for a build from a development
# source tree based on git repository information.
# It is assumed that by default the script is run in cmake script mode.
# If *not* called in script mode but used in generating cache variables,
# GEN_VERSION_INFO_INTERNAL has to be set ON.
#
# The following variables have to be previously defined:
# GIT_EXECUTABLE         - path to git binary (must be >=1.5.3)
# PROJECT_VERSION        - hard-coded version string (generated info is appended)
# PROJECT_SOURCE_DIR     - top level source directory (which has to be in git)
# VERSION_CMAKEIN        - path to an input template file
# VERSION_OUT            - path to the output file
# VERSION_NO_REMOTE_HASH - if set, GMX_GIT_REMOTE_HASH is not generated
#
# Output:
# i)  Script mode: VERSION_OUT is configured from the input VERSION_CMAKEIN
# using the variables listed below.
# ii) Cache variable mode: the variables below are set in cache.
#
# GMX_PROJECT_VERSION_STR   - version string
# GMX_GIT_HEAD_HASH         - git hash of current local HEAD
# GMX_GIT_REMOTE_HASH       - git hash of the first ancestor commit from the
#                             main Gromacs repository
#
# Szilard Pall (pszilard@cbr.su.se)
# Teemu Murtola (teemu.murtola@gmail.com)

if("${PROJECT_VERSION}" STREQUAL "")
    message(FATAL_ERROR "PROJECT_VERSION undefined!")
endif()

# if we're generating variables for cache unset the variables
if(GEN_VERSION_INFO_INTERNAL)
    set(GMX_PROJECT_VERSION_STR)
    set(GMX_GIT_HEAD_HASH)
    set(GMX_GIT_REMOTE_HASH)
endif()

# bail if the source tree is not in a git repository
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(FATAL_ERROR "Project source directory ${PROJECT_SOURCE_DIR} not in git")
endif()

# If git executable exists, build the development version string.
if(EXISTS "${GIT_EXECUTABLE}")
    # refresh git index
    execute_process(COMMAND ${GIT_EXECUTABLE} update-index -q --refresh
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        TIMEOUT 5
        OUTPUT_QUIET
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # get the full hash of the current HEAD
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE HEAD_HASH
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(GMX_GIT_HEAD_HASH ${HEAD_HASH})
    # extract the shortened hash (7 char)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE HEAD_HASH_SHORT
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # if there are local uncommitted changes, the build gets labeled "dirty"
    execute_process(COMMAND ${GIT_EXECUTABLE} diff-index --name-only HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE SRC_LOCAL_CHANGES
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT "${SRC_LOCAL_CHANGES}" STREQUAL "")
        set(DIRTY_STR "-dirty")
        set(GMX_GIT_HEAD_HASH "${GMX_GIT_HEAD_HASH} (dirty)")
    endif()

    # get the date of the HEAD commit
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-list -n1 "--pretty=format:%ci" HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE HEAD_DATE
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    string(REGEX REPLACE "\n| " ";" HEAD_DATE "${HEAD_DATE}")
    list(GET HEAD_DATE 2 HEAD_DATE)
    string(REGEX REPLACE "-" "" HEAD_DATE "${HEAD_DATE}")

    # compile the version string suffix
    set(VERSION_STR_SUFFIX "${HEAD_DATE}-${HEAD_HASH_SHORT}${DIRTY_STR}")

    if (DEFINED VERSION_NO_REMOTE_HASH)
        set(GMX_REMOTES "")
    else()
        # find the names of remotes that are located on the official gromacs
        # git/gerrit servers
        execute_process(COMMAND ${GIT_EXECUTABLE} config --get-regexp
                        "remote\\..*\\.url" "\\.gromacs\\.org[:/].*gromacs(\\.git)?$"
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GMX_REMOTES
            ERROR_VARIABLE EXEC_ERR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()

    # if there are remotes from the gromacs git servers, try to find ancestor
    # commits of the current HEAD from this remote;
    # otherwise, label the build "unknown"
    if("${GMX_REMOTES}" STREQUAL "")
        if (NOT DEFINED VERSION_NO_REMOTE_HASH)
            set(VERSION_STR_SUFFIX "${VERSION_STR_SUFFIX}-unknown")
        endif()
        set(GMX_GIT_REMOTE_HASH "unknown")
    else()
        string(REPLACE "\n" ";" GMX_REMOTES ${GMX_REMOTES})
        # construct a command pipeline that produces a reverse-time-ordered
        # list of commits and their annotated names in GMX_REMOTES
        # the max-count limit is there to put an upper bound on the execution time
        set(BASEREVCOMMAND "COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1000 HEAD")
        foreach(REMOTE ${GMX_REMOTES})
            string(REGEX REPLACE "remote\\.(.*)\\.url.*" "\\1" REMOTE ${REMOTE})
            set(BASEREVCOMMAND "${BASEREVCOMMAND} COMMAND ${GIT_EXECUTABLE} name-rev --stdin --refs=refs/remotes/${REMOTE}/*")
        endforeach(REMOTE)
        # this is necessary for CMake to properly split the variable into
        # parameters for execute_process().
        string(REPLACE " " ";" BASEREVCOMMAND ${BASEREVCOMMAND})
        # find the first ancestor in the list provided by rev-list (not
        # necessarily the last though) which is in GMX_REMOTE, extract the
        # hash and the number of commits HEAD is ahead with
        execute_process(${BASEREVCOMMAND}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE ANCESTOR_LIST
        )
        string(REGEX REPLACE "\n" ";" ANCESTOR_LIST "${ANCESTOR_LIST}")

        set(AHEAD 0)
        set(GMX_GIT_REMOTE_HASH "")
        foreach(ANCESTOR ${ANCESTOR_LIST})
            string(REPLACE "\n" "" HASH_AND_REVNAMES "${ANCESTOR}")
            string(REPLACE " " ";" HASH_AND_REVNAMES "${HASH_AND_REVNAMES}")
            list(LENGTH HASH_AND_REVNAMES COUNT)
            # stop and set the hash if we have a hit, otherwise loop and count
            # how far ahead is the local repo
            if(COUNT GREATER 1)
                LIST(GET HASH_AND_REVNAMES 0 GMX_GIT_REMOTE_HASH)
                break()
            endif()
            math(EXPR AHEAD ${AHEAD}+1)
        endforeach(ANCESTOR)
        # mark the build "local" if didn't find any commits that are from
        # remotes/${GMX_REMOTE}/*
        if("${GMX_GIT_REMOTE_HASH}" STREQUAL "")
            set(GMX_GIT_REMOTE_HASH "unknown")
            set(VERSION_STR_SUFFIX "${VERSION_STR_SUFFIX}-local")
        # don't print the remote hash if there are no local commits
        elseif("${GMX_GIT_REMOTE_HASH}" STREQUAL "${HEAD_HASH}")
            set(GMX_GIT_REMOTE_HASH "")
        else()
            set(GMX_GIT_REMOTE_HASH "${GMX_GIT_REMOTE_HASH} (${AHEAD} newer local commits)")
        endif()
    endif()

    # compile final version string
    set(GMX_PROJECT_VERSION_STR "${PROJECT_VERSION}-${VERSION_STR_SUFFIX}")
else()
    # the version has to be defined
    # set the GIT related information to "unknown"
    message(WARNING "Source tree seems to be a repository, but no compatible git is available, using hard-coded version string")
    set(GMX_PROJECT_VERSION_STR "${PROJECT_VERSION}")
    set(GMX_GIT_HEAD_HASH "unknown")
    set(GMX_GIT_REMOTE_HASH "unknown")
endif()

# if we're generating cache variables set these
# otherwise it's assumed that it's called in script mode to generate a file
if(GEN_VERSION_INFO_INTERNAL)
    set(GMX_PROJECT_VERSION_STR ${GMX_PROJECT_VERSION_STR}
        CACHE STRING "Gromacs version string" FORCE)
    set(GMX_GIT_HEAD_HASH ${GMX_GIT_HEAD_HASH}${DIRTY_STR}
        CACHE STRING "Current git HEAD commit object" FORCE)
    set(GMX_GIT_REMOTE_HASH ${GMX_GIT_REMOTE_HASH}
        CACHE STRING "Commmit object of the nearest ancestor present in the Gromacs git repository" FORCE)
    mark_as_advanced(GMX_GIT_HEAD_HASH GMX_GIT_REMOTE_HASH)
else()
    if("${VERSION_CMAKEIN}" STREQUAL "")
        message(FATAL_ERROR "Missing input parameter VERSION_CMAKEIN!")
    endif()
    if("${VERSION_OUT}" STREQUAL "")
        message(FATAL_ERROR "Missing input parameter VERSION_OUT!")
    endif()
    # Generate the output file.
    configure_file(${VERSION_CMAKEIN} ${VERSION_OUT})
endif()
