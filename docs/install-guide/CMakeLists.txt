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

# Note that the install-guide target, etc. can still be built
# independent of the webpage target.
if(MARKDOWN_CONFIGURE_IS_POSSIBLE AND PANDOC_EXECUTABLE)
    set(name "install-guide")
    configure_markdown(${name}.md)
    make_markdown_html(${name} ${CMAKE_CURRENT_BINARY_DIR}/${name}.md)
    make_markdown_pdf(${name} ${CMAKE_CURRENT_BINARY_DIR}/${name}.md)

    # Make the INSTALL file for CPack for the tarball. This gets put
    # into the tarball via the CPack rules in the top-level
    # CMakeLists.txt
    add_custom_command(
        OUTPUT final/INSTALL
        COMMAND ${CMAKE_COMMAND} -E make_directory final
        COMMAND ${PANDOC_EXECUTABLE} -t plain -o final/INSTALL ${name}.md
        DEPENDS
            ${CMAKE_CURRENT_BINARY_DIR}/${name}.md
        VERBATIM
        )

    # Add a top-level target for the webpage build to hook onto
    add_custom_target(${name}
        DEPENDS
            final/INSTALL
            ${HTML_OUTPUT_DIR}/${name}.html
            ${HTML_OUTPUT_DIR}/${name}.pdf
        VERBATIM
        )
endif()
