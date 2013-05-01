#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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

set(COULOMBINTEGRALS_SOURCES
     coulombintegrals/DSlater_1S_1S.cpp
     coulombintegrals/DSlater_1S_2S.cpp
     coulombintegrals/DSlater_1S_3S.cpp
     #coulombintegrals/DSlater_1S_4S.cpp
     #coulombintegrals/DSlater_1S_5S.cpp
     #coulombintegrals/DSlater_1S_6S.cpp
     coulombintegrals/DSlater_2S_2S.cpp
     coulombintegrals/DSlater_2S_3S.cpp
     #coulombintegrals/DSlater_2S_4S.cpp
     #coulombintegrals/DSlater_2S_5S.cpp
     #coulombintegrals/DSlater_2S_6S.cpp
     coulombintegrals/DSlater_3S_3S.cpp
     #coulombintegrals/DSlater_3S_4S.cpp
     #coulombintegrals/DSlater_3S_5S.cpp
     #coulombintegrals/DSlater_3S_6S.cpp
     #coulombintegrals/DSlater_4S_4S.cpp
     #coulombintegrals/DSlater_4S_5S.cpp
     #coulombintegrals/DSlater_4S_6S.cpp
     #coulombintegrals/DSlater_5S_5S.cpp
     #coulombintegrals/DSlater_5S_6S.cpp
     #coulombintegrals/DSlater_6S_6S.cpp
     coulombintegrals/Slater_1S_1S.cpp
     coulombintegrals/Slater_1S_2S.cpp
     coulombintegrals/Slater_1S_3S.cpp
     #coulombintegrals/Slater_1S_4S.cpp
     #coulombintegrals/Slater_1S_5S.cpp
     #coulombintegrals/Slater_1S_6S.cpp
     coulombintegrals/Slater_2S_2S.cpp
     coulombintegrals/Slater_2S_3S.cpp
     #coulombintegrals/Slater_2S_4S.cpp
     #coulombintegrals/Slater_2S_5S.cpp
     #coulombintegrals/Slater_2S_6S.cpp
     coulombintegrals/Slater_3S_3S.cpp
     #coulombintegrals/Slater_3S_4S.cpp
     #coulombintegrals/Slater_3S_5S.cpp
     #coulombintegrals/Slater_3S_6S.cpp
     #coulombintegrals/Slater_4S_4S.cpp
     #coulombintegrals/Slater_4S_5S.cpp
     #coulombintegrals/Slater_4S_6S.cpp
     #coulombintegrals/Slater_5S_5S.cpp
     #coulombintegrals/Slater_5S_6S.cpp
     #coulombintegrals/Slater_6S_6S.cpp
     coulombintegrals/gaussian_integrals.c
     coulombintegrals/slater_integrals.cpp
)

set(LIBGROMACS_SOURCES
    ${LIBGROMACS_SOURCES} ${COULOMBINTEGRALS_SOURCES} PARENT_SCOPE)

