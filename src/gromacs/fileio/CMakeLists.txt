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

file(GLOB FILEIO_SOURCES *.cpp *.c)

if(GMX_USE_PLUGINS)
  set(FILEIO_SOURCES ${FILEIO_SOURCES} ${CMAKE_SOURCE_DIR}/src/external/vmd_molfile/vmddlopen.c)
else()
  list(REMOVE_ITEM FILEIO_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/vmdio.c)
endif()

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${FILEIO_SOURCES} PARENT_SCOPE)

gmx_install_headers(
    confio.h
    enxio.h
    filenm.h
    gmxfio.h
    matio.h
    mdoutf.h
    mtxio.h
    pdbio.h
    tpxio.h
    trajectory_writing.h
    trnio.h
    trx.h
    trxio.h
    xdr_datatype.h
    xtcio.h
    xvgr.h
    )

if (GMX_USE_TNG AND TNG_IO_DEFINITIONS)
    set_property(SOURCE tngio.cpp
                 APPEND PROPERTY COMPILE_DEFINITIONS ${TNG_IO_DEFINITIONS})
    set_property(SOURCE tngio_for_tools.cpp
                 APPEND PROPERTY COMPILE_DEFINITIONS ${TNG_IO_DEFINITIONS})
endif()

if (BUILD_TESTING)
     add_subdirectory(tests)
endif()
