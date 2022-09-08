#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
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

gmx_add_unit_test(UtilityUnitTests utility-test
    CPP_SOURCE_FILES
        alignedallocator.cpp
        arrayref.cpp
        booltype.cpp
        bitmask32.cpp bitmask64.cpp bitmask128.cpp
        cstringutil.cpp
        defaultinitializationallocator.cpp
        enumerationhelpers.cpp
        fixedcapacityvector.cpp
        inmemoryserializer.cpp
        keyvaluetreeserializer.cpp
        keyvaluetreetransform.cpp
        listoflists.cpp
        logger.cpp
        message_string_collector.cpp
        path.cpp
        physicalnodecommunicator.cpp
        range.cpp
        strconvert.cpp
        stringtoenumvalueconverter.cpp
        stringutil.cpp
        template_mp.cpp
        textreader.cpp
        textwriter.cpp
        typetraits.cpp
        )
# TODO: Remove `legacy_modules` once specific modules are explicitly linked.
target_link_libraries(utility-test PRIVATE legacy_modules)
target_link_libraries(utility-test PRIVATE utility)

gmx_add_mpi_unit_test(UtilityMpiUnitTests utility-mpi-test 2
    CPP_SOURCE_FILES
        physicalnodecommunicator_mpi.cpp
        )
