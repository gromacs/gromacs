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

# Propagate all flags passed to parent find_package() to the config call below.
set(_gmx_find_args "")
if (GROMACS_FIND_VERSION)
    list(APPEND _gmx_find_args ${GROMACS_FIND_VERSION})
    if (GROMACS_FIND_VERSION_EXACT)
        list(APPEND _gmx_find_args EXACT)
    endif()
endif()
if (GROMACS_FIND_REQUIRED)
    list(APPEND _gmx_find_args REQUIRED)
endif()
if (GROMACS_FIND_QUIETLY)
    list(APPEND _gmx_find_args QUIET)
endif()
# Delegate all the actual work to the package configuration files, except
# that select the configuration file based on suffix.
# TODO: Check that these produce sensible error messages and improve if possible.
if (DEFINED GROMACS_SUFFIX)
    find_package(GROMACS ${_gmx_find_args} CONFIG NAMES gromacs${GROMACS_SUFFIX})
else()
    find_package(GROMACS ${_gmx_find_args} CONFIG)
endif()
unset(_gmx_find_args)
