#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

# Managing configuration for Intel Threading Building Blocks

if(GMX_TBB)

    if(TBB_FOUND)
        include_directories(${TBB_INCLUDE_DIRS})

        # Base TBB library is available if TBB_FOUND=TRUE, but tbbmalloc is optional
        list(APPEND GMX_EXTRA_LIBRARIES ${TBB_LIBRARIES})
        if(TBB_MALLOC_LIBRARIES)
            list(APPEND GMX_EXTRA_LIBRARIES ${TBB_MALLOC_LIBRARIES})
        endif()

        if(APPLE)
            # The commercial version of TBB is screwed up on OS X, and does not
            # set its rpath id correctly. We fix this by setting rpath explicitly
            # here, and then alter the build targets to add the rpath prefix.
            # The last element of list should contain a filename with full path
            list(GET TBB_LIBRARIES -1 TMP_TBB_PATH)
            # get the directory part and add to rpath
            get_filename_component(TMP_TBB_PATH "${TMP_TBB_PATH}" DIRECTORY)
            list(APPEND CMAKE_INSTALL_RPATH ${TMP_TBB_PATH})
            list(APPEND CMAKE_BUILD_RPATH ${TMP_TBB_PATH})
        endif()
    else()
        message(FATAL_ERROR "Intel TBB support requested, but not found.")
    endif()
endif()
