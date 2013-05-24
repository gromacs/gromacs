#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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
    set(${NAME} ${DEFAULT} CACHE STRING ${_description})
    set_property(CACHE ${NAME} PROPERTY STRINGS ${_allowed})

    # Check that the value is one of the allowed
    set(_org_value "${${NAME}}")
    string(TOUPPER "${${NAME}}" ${NAME})
    string(TOUPPER "${_allowed}" _allowed_as_upper)
    list(FIND _allowed_as_upper "${${NAME}}" _found_index)
    if (_found_index EQUAL -1)
        message(FATAL_ERROR "Invalid value for ${NAME}: ${_org_value}.  "
                            "Pick one of: ${_allowed_comma_separated}")
    endif ()
    # Always provide the upper-case value to the caller
    set(${NAME} "${${NAME}}" PARENT_SCOPE)
endfunction()

# Convenience function for reporting a fatal error for an invalid input value
function(GMX_INVALID_OPTION_VALUE NAME)
    message(FATAL_ERROR "Invalid value for ${NAME}: ${${NAME}}")
endfunction()
