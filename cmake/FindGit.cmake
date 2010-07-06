# Try to find the git version control tool
# 
# The module defines the following variables:
# 
# Git_EXECUTABLE    - path the the git executable
# Git_VERSION       - git version
# Git_FOUND         - tru if git was found, false otherwise
#
# Author: Szilard Pall (pszilard@cbr.su.se)

if(Git_EXECUTABLE AND Git_VERSION)
    set(Git_FIND_QUIETLY TRUE)
endif()

# search for git binary
find_program(Git_EXECUTABLE git
    PATHS ENV PATH
    CACHE DOC "Git version control tool")

if(NOT Git_EXECUTABLE)
    set(_err_msg "Git executable not found")
    if(Git_FIND_REQUIRED)
        message(FATAL_ERROR " ${_err_msg}")
    elseif(NOT Git_FIND_QUIETLY)
        message("${_err_msg}")
    endif()
endif()

# parse version
if(Git_EXECUTABLE AND NOT Git_VERSION)
    execute_process(COMMAND ${Git_EXECUTABLE} "--version"
        OUTPUT_VARIABLE _exec_out
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "git version (.*)" "\\1" Git_VERSION ${_exec_out})
    set(Git_VERSION ${Git_VERSION} CACHE STRING "Git version")
    
    # check version
    set(_git_version_ok TRUE)
    # this should at some point become VERSION_EQUAL
    if(Git_FIND_VERSION_EXACT AND NOT Git_VERSION STREQUAL Git_FIND_VERSION)
        set(_err_msg "Found git version ${Git_VERSION} but this does not match the requested ${Git_FIND_VERSION}")
        if(Git_FIND_REQUIRED)
            message(FATAL_ERROR " ${_err_msg}")
        elseif(NOT Git_FIND_QUIETLY)
            message("${_err_msg}")
        endif()
        set(_git_version_ok FALSE)
    endif()
    # this should at some point become VERSION_LESS
    if(Git_FIND_VERSION AND Git_VERSION STRLESS Git_FIND_VERSION)
        set(_err_msg "Found git version ${Git_VERSION} but this is less then the requested ${Git_FIND_VERSION}")
        if(Git_FIND_REQUIRED)
            message(FATAL_ERROR " ${_err_msg}")
        elseif(NOT Git_FIND_QUIETLY)
            message("${_err_msg}")
        endif()
        set(_git_version_ok FALSE)
    endif()

endif()
set(Git_FOUND True)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Git DEFAULT_MSG 
    Git_EXECUTABLE 
    Git_VERSION
    _git_version_ok)

mark_as_advanced(Git_EXECUTABLE Git_VERSION)
