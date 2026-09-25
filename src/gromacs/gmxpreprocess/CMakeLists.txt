#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2011- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# Set up the module library
add_library(gmxpreprocess INTERFACE)

set(GMXPREPROCESS_SOURCES
    add_par.cpp
    calch.cpp
    compute_io.cpp
    convparm.cpp
    editconf.cpp
    fflibutil.cpp
    gen_ad.cpp
    gen_maxwell_velocities.cpp
    gen_vsite.cpp
    genconf.cpp
    genhydro.cpp
    genion.cpp
    genrestr.cpp
    gmxcpp.cpp
    gpp_atomtype.cpp
    gpp_bond_atomtype.cpp
    gpp_nextnb.cpp
    grompp.cpp
    h_db.cpp
    hackblock.cpp
    hizzie.cpp
    insert_molecules.cpp
    makeexclusiondistances.cpp
    massrepartitioning.cpp
    nm2type.cpp
    pdb2gmx.cpp
    pdb2top.cpp
    pgutil.cpp
    readir.cpp
    readpull.cpp
    readrot.cpp
    resall.cpp
    solvate.cpp
    specbond.cpp
    ter_db.cpp
    tomorse.cpp
    topdirs.cpp
    topio.cpp
    toppush.cpp
    topshake.cpp
    toputil.cpp
    vsite_parm.cpp
    x2top.cpp
    xlate.cpp
    )
list(TRANSFORM GMXPREPROCESS_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${GMXPREPROCESS_SOURCES} PARENT_SCOPE)

# Source files have the following private module dependencies.
target_link_libraries(gmxpreprocess PRIVATE
                      #                      gmxlib
                      #                      math
                      #                      mdtypes
                      #                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(gmxpreprocess PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(gmxpreprocess PUBLIC
target_link_libraries(gmxpreprocess INTERFACE
                      legacy_api
                      mdtypes
                      )

# TODO: when is an OBJECT target
#target_link_libraries(gmxpreprocess PUBLIC legacy_api)
#target_link_libraries(gmxpreprocess PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(gmxpreprocess PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(gmxpreprocess PRIVATE legacy_modules)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
