# Generate Gromacs development build version information.
#
# The script generates version information for a build from a development 
# source tree based on git repository information. 
# It is assumed that by default the script is run in cmake script mode.
# If *not* called in script mode but used in generating cache variables,  
# GEN_VERSION_INFO_INTERNAL has to be set ON.
#
# The following variables have to be previously defined: 
# PROJECT_VERSION       - hard-coded version string, should have the following structure: 
#                       VERSION[-dev-SUFFIX] where the VERSION can have any form and the suffix 
#                       is optional but should start with -dev
# PROJECT_SOURCE_DIR    - top level source directory 
# VERSION_C_CMAKEIN     - path to the version.c.cmakein file 
# VERSION_C_OUT         - path to the version.c output file
#
# Output: 
# i)  Script mode: version.c configured from the input version.c.cmakein using 
# the variables listed below. 
# ii) Cache variable mode: the varables below are set in cache.
#
# GMX_PROJECT_VERSION_STR   - version string 
# GMX_GIT_HEAD_HASH         - git hash of current local HEAD 
# GMX_GIT_REMOTE_HASH       - git hash of the first ancestor commit from the 
#                             main Gromacs repository 
#

if(${PROJECT_VERSION} STREQUAL "")
    message(FATAL_ERROR "PROJECT_VERSION undefined!")
endif()
set(VER ${PROJECT_VERSION})

# if we're generating variables for cache unset the variables 
if(GEN_VERSION_INFO_INTERNAL)
    set(GMX_PROJECT_VERSION_STR)
    set(GMX_GIT_HEAD_HASH)
    set(GMX_GIT_REMOTE_HASH)
endif()

set(GIT_BIN)
find_program(GIT_BIN "git")
mark_as_advanced(GIT_BIN)

# if the source tree is a git repository extract the necessary information for  
# building the development version string 
if(NOT "${GIT_BIN}" MATCHES ".*-NOTFOUND" 
   AND IS_DIRECTORY "${PROJECT_SOURCE_DIR}/.git")
    # refresh git index 
    execute_process(COMMAND ${GIT_BIN} update-index -q --refresh
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        TIMEOUT 5
        OUTPUT_QUIET
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    ) 

   # get the full hash of the current HEAD 
    execute_process(COMMAND ${GIT_BIN} rev-parse -q HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        TIMEOUT 5
        OUTPUT_VARIABLE GMX_GIT_HEAD_HASH
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # extract the shortened hash (7 char)
    string(SUBSTRING ${GMX_GIT_HEAD_HASH} 0 5 HEAD_HASH_SHORT) 

    # if there are local uncommitted changes, the build gets labeled "dirty"
    execute_process(COMMAND ${GIT_BIN} diff-index --name-only HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        TIMEOUT 5
        OUTPUT_VARIABLE SRC_LOCAL_CHANGES
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )   
    if(NOT ${SRC_LOCAL_CHANGES} STREQUAL "")
        set(DIRTY_STR "-dirty")
        set(GMX_GIT_HEAD_HASH "${GMX_GIT_HEAD_HASH} (dirty)")
    endif()

    # get the date of the HEAD commit
    execute_process(COMMAND ${GIT_BIN} rev-list -n1 "--pretty=format:%ci" HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        TIMEOUT 5
        OUTPUT_VARIABLE HEAD_DATE
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    string(REGEX REPLACE "\n| " ";" HEAD_DATE ${HEAD_DATE})
    list(GET HEAD_DATE 2 HEAD_DATE)
    string(REGEX REPLACE "-" "" HEAD_DATE ${HEAD_DATE})

    # compile the version string suffix
    set(VERSION_STR_SUFFIX "${HEAD_DATE}-${HEAD_HASH_SHORT}${DIRTY_STR}") 
    
    # find the name of the remote which is located on the official gromacs git server
    execute_process(COMMAND ${GIT_BIN} config --get-regexp 
                    "remote\\..*\\.url" "git\\.gromacs\\.org[:|/]gromacs"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        TIMEOUT 5
        OUTPUT_VARIABLE GMX_REMOTE
        ERROR_VARIABLE EXEC_ERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # if there's a remote from the gromacs git, try to find ancestor commits of the 
    # current HEAD from this remote; otherwise, label the buld "unknown"
    if(GMX_REMOTE STREQUAL "")
        set(VERSION_STR_SUFFIX "${VERSION_STR_SUFFIX}-unknown")
        set(GMX_GIT_REMOTE_HASH "unknown")        
    else()         
        string(REGEX REPLACE "remote\\.(.*)\\.url.*" "\\1" GMX_REMOTE ${GMX_REMOTE})
        # find the first ancestor in the list provided by rev-list (not 
        # necessarily the last though) which is in GMX_REMOTE, extract the 
        # hash and the number of commits HEAD is ahead with 
        execute_process(COMMAND ${GIT_BIN} rev-list --max-count=100 HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            TIMEOUT 5
            OUTPUT_VARIABLE ANCESTOR_LIST
        )
        string(REGEX REPLACE "\n" ";" ANCESTOR_LIST ${ANCESTOR_LIST})

        set(AHEAD 0)
        set(GMX_GIT_REMOTE_HASH "")
        foreach(OBJ ${ANCESTOR_LIST})
            execute_process(COMMAND ${GIT_BIN} name-rev --refs=refs/remotes/${GMX_REMOTE}/* ${OBJ}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                TIMEOUT 5
                OUTPUT_VARIABLE HASH_AND_REVNAME
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            string(REGEX REPLACE "\n" "" HASH_AND_REVNAME ${HASH_AND_REVNAME})
            string(REGEX REPLACE " " ";" HASH_AND_REVNAME ${HASH_AND_REVNAME})
            list(GET HASH_AND_REVNAME 0 GMX_GIT_REMOTE_HASH) 
            list(GET HASH_AND_REVNAME 1 REVNAME)
            # stop and set the hash if we have a hit, otherwise loop and count
            # how far ahead is the local repo
            if(${REVNAME} MATCHES "remotes/${GMX_REMOTE}/.*")
                set(GMX_GIT_REMOTE_HASH 
                        "${GMX_GIT_REMOTE_HASH} (${AHEAD} newer local commits)")
                break()
            else()
                math(EXPR AHEAD ${AHEAD}+1)
            endif()
        endforeach(OBJ)
        # mark the build "local" if didn't find any commits that are from 
        # remotes/${GMX_REMOTE}/*
        if(${GMX_GIT_REMOTE_HASH} STREQUAL "")
            set(GMX_GIT_REMOTE_HASH "unknown")
            set(VERSION_STR_SUFFIX "${VERSION_STR_SUFFIX}-local") 
        endif()
    endif()

    # compile final version string, if there is already a -dev suffix in VER 
    # remove everything after this and replace it with the generated suffix
    string(REGEX REPLACE "(.*)-dev.*" "\\1" VER ${VER})
    set(GMX_PROJECT_VERSION_STR "${VER}-dev-${VERSION_STR_SUFFIX}")
else()
    # the version has to be defined - if not we're not using version.h/.c and set 
    # the GIT related information to "unknown"
    message(WARNING " Not a git repository and/or git not found, using hard-coded version.")
    set(GMX_PROJECT_VERSION_STR "${PROJECT_VERSION}")
    set(GMX_GIT_HEAD_HASH "unknown")
    set(GMX_GIT_REMOTE_HASH "unknown")
    set(USE_VERSION_H OFF)
endif()

# if we're generating cache variables set these
# otherwise it's assumed that it's called in script mode to generate version.c 
if(GEN_VERSION_INFO_INTERNAL)
    set(GMX_PROJECT_VERSION_STR ${GMX_PROJECT_VERSION_STR} 
        CACHE STRING "Gromacs version string" FORCE)
    set(GMX_GIT_HEAD_HASH ${GMX_GIT_HEAD_HASH}${DIRTY_STR}  
        CACHE STRING "Current git HEAD commit object" FORCE)
    set(GMX_GIT_REMOTE_HASH ${GMX_GIT_REMOTE_HASH} 
        CACHE STRING "Commmit object of the nearest ancestor present in the Gromacs git repository" FORCE)
    mark_as_advanced(GMX_GIT_HEAD_HASH GMX_GIT_REMOTE_HASH)
else()
    if(${VERSION_C_CMAKEIN} STREQUAL "")
        message(FATAL_ERROR "Missing input parameter VERSION_C_CMAKEIN!")
    endif()
    if(${VERSION_C_OUT} STREQUAL "")
        message(FATAL_ERROR "Missing input parameter VERSION_C_OUT!")
    endif()
    # generate version.c
   configure_file(${VERSION_C_CMAKEIN} ${VERSION_C_OUT})    
endif()
