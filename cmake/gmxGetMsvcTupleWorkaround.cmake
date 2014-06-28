#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

# GMock uses tuples extensively, and MSVC bundles a tuple library that
# is not compatible with the standard. r675 of googletest works around
# this properly, but that's not in GMock 1.7.0. That logic is
# duplicated here. See
# https://code.google.com/p/googletest/source/detail?r=675#, but note
# that its summary does not represent its code correctly.
#
# This function should be called to get the compile definitions
# suitable for working around MSVC to compile GMock, if any.
# Returns a string of options in VARIABLE
function(GET_MSVC_TUPLE_WORKAROUND_DEFINITIONS VARIABLE)
    set(${VARIABLE} "")
    if(MSVC)
        if(MSVC_VERSION VERSION_LESS 1600)
            list(APPEND ${VARIABLE} "/D GTEST_USE_OWN_TR1_TUPLE=1")
        else()
            list(APPEND ${VARIABLE} "/D GTEST_USE_OWN_TR1_TUPLE=0")
            if(MSVC_VERSION VERSION_EQUAL 1700)
                # Fixes Visual Studio 2012
                list(APPEND ${VARIABLE} "/D _VARIADIC_MAX=10")
            endif()
        endif()
    endif()
    set(${VARIABLE} ${${VARIABLE}} PARENT_SCOPE)
endfunction()
