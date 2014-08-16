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

# GMock uses tuples extensively, but since it is part of tr1 we cannot
# assume it is present - the Fujitsu compilers on K computer is one
# example of a system where it is not, and this is likely the case for many
# embedded systems too.
#
# In general, we can ask gmock to use its own internal implementation by
# defining GTEST_USE_OWN_TR1_TUPLE=1. This is safe for Gromacs, since we do
# not use tr1/tuple.h elsewhere. However, this workaround itself will not
# compile on MSVC - but this is the only architecture we know where it does
# not work. To make things even worse, on MSVC even the standard tr1/tuple
# implementation is not fully compatible, but we can make it work by setting
# _VARIADIC_MAX=10. This is similar to r675 of googletest, which is still not
# present in GMock 1.7.0. See
# https://code.google.com/p/googletest/source/detail?r=675#, but note
# that its summary does not represent its code correctly.
#
# To make all this work without user intervention, we first check what compiler
# we have, and if it is not MSVC we simply use the internal version. If we are
# using MSVC, we rely on the compiler's version, but set the variable necessary
# to make it compatible.
#
# This function should be called to get the compile definitions
# suitable for making sure we have a tuple implementation that works with gmock.
#
# Returns a string of options in VARIABLE
function(GET_GMOCK_TUPLE_WORKAROUND VARIABLE)
    set(${VARIABLE} "")
    if(MSVC OR (WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "Intel"))
        if(MSVC_VERSION VERSION_EQUAL 1700)
            # Fixes Visual Studio 2012
            set(${VARIABLE} "_VARIADIC_MAX=10")
        endif()
        # For MSVC different from 2012, or Intel on Windows, we reply on tr1/tuple working.
    else()
        # Use the built-in version on all other compilers.
        set(${VARIABLE} "GTEST_USE_OWN_TR1_TUPLE=1")
    endif()
    set(${VARIABLE} ${${VARIABLE}} PARENT_SCOPE)
endfunction()
