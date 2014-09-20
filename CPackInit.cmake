#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#TODO: add check that source doesn't contain any untracked files
if(NOT CPACK_INSTALL_CMAKE_PROJECTS) #building source package
    get_filename_component(CMAKE_BINARY_DIR ${CPACK_OUTPUT_CONFIG_FILE} PATH)
    if (NOT EXISTS "${CMAKE_BINARY_DIR}/docs/man/man1/gmx-view.1" OR
        NOT EXISTS "${CMAKE_BINARY_DIR}/docs/install-guide/final/INSTALL" OR
        NOT EXISTS "${CMAKE_BINARY_DIR}/docs/old-html/final/online.html" OR
        NOT EXISTS "${CMAKE_BINARY_DIR}/src/programs/completion/gmx-completion.bash" OR
        NOT EXISTS "${CMAKE_BINARY_DIR}/docs/old-html/final/programs/gmx-view.html")
        message(FATAL_ERROR
            "To create a complete source package, bash completions, "
            "man and HTML pages, and INSTALL need to be generated. "
            "Run 'make completion man html install-guide' to build "
            "these parts. You can also configure with "
            "GMX_BUILD_HELP=ON to automatically build the HTML parts.")
    endif()
else()
    if (NOT CPACK_GMX_BUILD_HELP)
        message(WARNING
            "To create a complete binary package, bash completions, and "
            "man and HTML pages need to be generated. "
            "You need to configure with GMX_BUILD_HELP=ON to include all "
            "in the binary package.")
        # Building the man, html, ... targets is not sufficient because than the
        # install is still not done.
    endif()
endif()
