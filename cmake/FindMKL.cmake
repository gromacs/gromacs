#
# This file is part of the GROMACS molecular simulation package,
# version 4.6
#
# Copyright (c) 2012, by the GROMACS development team, led by
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
# - Find Intel MKL
# Find the Intel Math Kernel Library, version 6.0 or later
#
#  MKL_INCLUDE_DIR - where to find mkl_dfti.h
#  MKL_LIBRARIES   - List of libraries when using MKL
#  MKL_FOUND       - True if MKL was found

if (MKL_INCLUDE_DIR)
  # Already in cache, be silent
  set (MKL_FIND_QUIETLY TRUE)
endif (MKL_INCLUDE_DIR)

find_path (MKL_INCLUDE_DIR mkl_dfti.h)
find_library (MKL_LIBRARIES mkl_core)

# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIR)

# MKL Libraries change name ALL the time, so the user will have to set these...
# mark_as_advanced (MKL_LIBRARIES MKL_INCLUDE_DIR)


