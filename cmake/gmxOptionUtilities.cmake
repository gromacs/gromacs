#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014, by the GROMACS development team, led by
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

# Helper functions for managing more complex options
#

# Creates a string cache variable with multiple choices
#
# Usage:
#   gmx_option_multichoice(VAR "Description" DEFAULT_VALUE
#                          Value1 Value2 ... ValueN)
# Output:
#   VAR is set in the cache and in the caller's scope. The caller can assume
#   that it is always one of the provided values, converted to uppercase.
#
# Main benefit is that the list of allowed values only needs to be provided
# once, and gets used in multiple contexts:
#   1. It is automatically added to the description.
#   2. It is set as the STRINGS property of the created cache variable for use
#      with CMake GUIs.
#   3. The user-provided value is checked against the list, and a fatal error
#      is produced if the value is not known.  The caller does not need to
#      produce good error messages in cases where it may be necessary to check
#      for the validity again.
# As a special case, any "[built-in]" string in the allowed values is ignored
# when checking the user-provided value, but is added to all user-visible
# messages.
#
# It appears that ccmake does not use the STRINGS property, but perhaps some
# day...
#
function(GMX_OPTION_MULTICHOICE NAME DESCRIPTION DEFAULT)
    # Some processing of the input values
    string(REPLACE ";" ", " _allowed_comma_separated "${ARGN}")
    set(_description "${DESCRIPTION}. Pick one of: ${_allowed_comma_separated}")
    string(REPLACE "[built-in]" "" _allowed "${ARGN}")

    # Set the cache properties
    set(${NAME} ${DEFAULT} CACHE STRING "${_description}")
    set_property(CACHE ${NAME} PROPERTY STRINGS ${_allowed})

    # Check that the value is one of the allowed
    set(_org_value "${${NAME}}")
    string(TOUPPER "${${NAME}}" ${NAME})
    string(TOUPPER "${_allowed}" _allowed_as_upper)
    list(FIND _allowed_as_upper "${${NAME}}" _found_index)
    if (_found_index EQUAL -1)
        message(FATAL_ERROR "Invalid value for ${NAME}: ${_org_value}.  "
                            "Pick one of: ${_allowed_comma_separated}")
    endif()
    # Always provide the upper-case value to the caller
    set(${NAME} "${${NAME}}" PARENT_SCOPE)
endfunction()

# Convenience function for reporting a fatal error for an invalid input value
function(GMX_INVALID_OPTION_VALUE NAME)
    message(FATAL_ERROR "Invalid value for ${NAME}: ${${NAME}}")
endfunction()

# Declares a cache variable with ON/OFF/AUTO values
#
# Usage:
#   gmx_option_trivalue(VAR "Description" DEFAULT)
#
# Output:
#   VAR is created in the cache, and the caller can assume that the value is
#   always one of ON/OFF/AUTO.  Additionally, VAR_AUTO is set if value is AUTO,
#   and VAR_FORCE is set if value is ON.
#   These make it convenient to check for any combination of states with simple
#   if() statements (simple if(VAR) matches AUTO and ON).
function(GMX_OPTION_TRIVALUE NAME DESCRIPTION DEFAULT)
    set(_description "${DESCRIPTION}. ON/OFF/AUTO")
    set(${NAME} ${DEFAULT} CACHE STRING "${_description}")
    set_property(CACHE ${NAME} PROPERTY STRINGS ON OFF AUTO)

    set(${NAME}_AUTO OFF)
    set(${NAME}_FORCE OFF)
    string(TOUPPER "${${NAME}}" ${NAME})
    if ("${${NAME}}" STREQUAL "AUTO")
        set(${NAME}_AUTO ON)
    elseif (${NAME})
        set(${NAME}_FORCE ON)
        set(${NAME} ON)
    else()
        set(${NAME} OFF)
    endif()
    # Always provide the sanitized value to the caller
    set(${NAME}       "${${NAME}}"       PARENT_SCOPE)
    set(${NAME}_AUTO  "${${NAME}_AUTO}"  PARENT_SCOPE)
    set(${NAME}_FORCE "${${NAME}_FORCE}" PARENT_SCOPE)
endfunction()

# Hides or shows a cache value based on conditions
#
# Usage:
#   gmx_add_cache_dependency(VAR TYPE CONDITIONS VALUE)
# where
#   VAR        is a name of a cached variable
#   TYPE       is the type of VAR
#   CONDITIONS is a list of conditional expressions (see below)
#   VALUE      is a value that is set to VAR if CONDITIONS is not satisfied
#
# Evaluates each condition in CONDITIONS, and if any of them is false,
# VAR is marked internal (hiding it from the user) and its value is set to
# VALUE.  The previous user-set value of VAR is still remembered in the cache,
# and used when CONDITIONS become true again.
#
# The conditions is a semicolon-separated list of conditions as specified for
# CMake if() statements, such as "GMX_FFT_LIBRARY STREQUAL FFTW3",
# "NOT GMX_MPI" or "GMX_MPI;NOT GMX_DOUBLE".  Note that quotes within the
# expressions don't work for some reason (even if escaped).
#
# The logic is adapted from cmake_dependent_option().
#
function(GMX_ADD_CACHE_DEPENDENCY NAME TYPE CONDITIONS FORCED_VALUE)
    set(_available TRUE)
    foreach(_cond ${CONDITIONS})
        string(REGEX REPLACE " +" ";" _cond_as_list ${_cond})
        if (${_cond_as_list})
        else()
            set(_available FALSE)
        endif()
    endforeach()
    if (_available)
        set_property(CACHE ${NAME} PROPERTY TYPE ${TYPE})
    else()
        set(${NAME} "${FORCED_VALUE}" PARENT_SCOPE)
        set_property(CACHE ${NAME} PROPERTY TYPE INTERNAL)
    endif()
endfunction()

# Works like cmake_dependent_option(), but allows for an arbitrary cache value
# instead of only an ON/OFF option
#
# Usage:
#   gmx_dependent_cache_variable(VAR "Description" TYPE DEFAULT CONDITIONS)
#
# Creates a cache variable VAR with the given description, type and default
# value.  If any of the conditions listed in CONDITIONS is not true, then
# the cache variable is marked internal (hiding it from the user) and the
# value of VAR is set to DEFAULT.  The previous user-set value of VAR is still
# remembered in the cache, and used when CONDITIONS become true again.
# Any further changes to the variable can be done with simple set()
# (or set_property(CACHE VAR PROPERTY VALUE ...) if the cache needs to be
# modified).
#
# See gmx_add_cache_dependency() on how to specify the conditions.
#
macro(GMX_DEPENDENT_CACHE_VARIABLE NAME DESCRIPTION TYPE DEFAULT CONDITIONS)
    set(${NAME} "${DEFAULT}" CACHE ${TYPE} "${DESCRIPTION}")
    gmx_add_cache_dependency(${NAME} ${TYPE} "${CONDITIONS}" "${DEFAULT}")
endmacro()

# Works like cmake_dependent_option(), but reuses the code from
# gmx_dependent_cache_variable() to make sure both behave the same way.
macro(GMX_DEPENDENT_OPTION NAME DESCRIPTION DEFAULT CONDITIONS)
    gmx_dependent_cache_variable(${NAME} "${DESCRIPTION}" BOOL "${DEFAULT}" "${CONDITIONS}")
endmacro()

# Checks if one or more cache variables have changed
#
# Usage:
#   gmx_check_if_changed(RESULT VAR1 VAR2 ... VARN)
#
# Sets RESULT to true if any of the given cache variables VAR1 ... VARN has
# changes since the last call to this function for that variable.
# Changes are tracked also across CMake runs.
function(GMX_CHECK_IF_CHANGED RESULT)
    set(_result FALSE)
    foreach (_var ${ARGN})
        if (NOT "${${_var}}" STREQUAL "${${_var}_PREVIOUS_VALUE}")
            set(_result TRUE)
        endif()
        set(${_var}_PREVIOUS_VALUE "${${_var}}" CACHE INTERNAL
            "Previous value of ${_var} for change tracking")
    endforeach()
    set(${RESULT} ${_result} PARENT_SCOPE)
endfunction()
