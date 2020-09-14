#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

gmx_add_libgromacs_sources(
    legacymdrunoptions.cpp
    legacysimulator.cpp
    md.cpp
    mdmodules.cpp
    membedholder.cpp
    mimic.cpp
    minimize.cpp
    replicaexchange.cpp
    rerun.cpp
    runner.cpp
    shellfc.cpp
    simulationcontext.cpp
    simulationinput.cpp
    simulationinputhandle.cpp
    simulatorbuilder.cpp
    tpi.cpp
    )

# TODO: Find a home for this header and a scheme for installation.
# This header straddles the installed libraries and is a transitive interface
# from libgromacs to libgmxapi to libgmxapi clients. Near term efforts are
# planned to consolidate these libraries and refine the public interface.
# In the mean time, this is kept out of the gmxapi to avoid circular dependencies
# between the libraries. It can be moved to a separate interface-only target,
# but that is beyond the scope of issue #3379, which introduces the header.
# Ref #3152 and #3652
install(FILES simulationinputhandle.h
        DESTINATION include/gromacs/mdrun)

if (BUILD_TESTING)
# TODO import this from src/programs/mdrun/tests
#    add_subdirectory(tests)
endif()
