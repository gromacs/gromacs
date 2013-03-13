#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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
if(X11_FOUND)

set(NGMX_COMMON_SOURCE
    xutil.c         xmb.c           x11.c
    xdlghi.c        fgrid.c         xdlgitem.c
    popup.c         pulldown.c      xdlg.c
    manager.c       buttons.c       nmol.c
    nleg.c          dialogs.c       logo.c
    filter.c        molps.c)

include_directories(${X11_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List of programs with single corresponding *.c source file,
# used to create build rules automatically.
#
set(NGMX_PROGRAMS 
    ngmx g_xrama)

foreach(PROG ${NGMX_PROGRAMS})
        add_executable(${PROG} ${PROG}.c ${NGMX_COMMON_SOURCE})
        gmx_add_man_page(${PROG})
        target_link_libraries(${PROG} gmx ${X11_LIBRARIES})
        set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG) 

install(TARGETS ${NGMX_PROGRAMS}
        COMPONENT ngmx
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(X11_FOUND)
