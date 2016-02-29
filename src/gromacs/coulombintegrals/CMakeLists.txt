#
# This file is part of the alexandria code base
#
# Copyright (C) 2014-2016 David van der Spoel and Paul J. van Maaren
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

set(COULOMB_SOURCES
     coulombintegrals/DSlater_1S_1S.cpp
     coulombintegrals/DSlater_1S_2S.cpp
     coulombintegrals/DSlater_1S_3S.cpp
     coulombintegrals/DSlater_1S_4S.cpp
     #coulombintegrals/DSlater_1S_5S.cpp
     #coulombintegrals/DSlater_1S_6S.cpp
     coulombintegrals/DSlater_2S_2S.cpp
     coulombintegrals/DSlater_2S_3S.cpp
     coulombintegrals/DSlater_2S_4S.cpp
     #coulombintegrals/DSlater_2S_5S.cpp
     #coulombintegrals/DSlater_2S_6S.cpp
     coulombintegrals/DSlater_3S_3S.cpp
     coulombintegrals/DSlater_3S_4S.cpp
     #coulombintegrals/DSlater_3S_5S.cpp
     #coulombintegrals/DSlater_3S_6S.cpp
     coulombintegrals/DSlater_4S_4S.cpp
     #coulombintegrals/DSlater_4S_5S.cpp
     #coulombintegrals/DSlater_4S_6S.cpp
     #coulombintegrals/DSlater_5S_5S.cpp
     #coulombintegrals/DSlater_5S_6S.cpp
     #coulombintegrals/DSlater_6S_6S.cpp
     coulombintegrals/Slater_1S_1S.cpp
     coulombintegrals/Slater_1S_2S.cpp
     coulombintegrals/Slater_1S_3S.cpp
     coulombintegrals/Slater_1S_4S.cpp
     #coulombintegrals/Slater_1S_5S.cpp
     #coulombintegrals/Slater_1S_6S.cpp
     coulombintegrals/Slater_2S_2S.cpp
     coulombintegrals/Slater_2S_3S.cpp
     coulombintegrals/Slater_2S_4S.cpp
     #coulombintegrals/Slater_2S_5S.cpp
     #coulombintegrals/Slater_2S_6S.cpp
     coulombintegrals/Slater_3S_3S.cpp
     coulombintegrals/Slater_3S_4S.cpp
     #coulombintegrals/Slater_3S_5S.cpp
     #coulombintegrals/Slater_3S_6S.cpp
     coulombintegrals/Slater_4S_4S.cpp
     #coulombintegrals/Slater_4S_5S.cpp
     #coulombintegrals/Slater_4S_6S.cpp
     #coulombintegrals/Slater_5S_5S.cpp
     #coulombintegrals/Slater_5S_6S.cpp
     #coulombintegrals/Slater_6S_6S.cpp
     coulombintegrals/gaussian_integrals.c
     coulombintegrals/slater_integrals.cpp
	)

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${COULOMB_SOURCES} PARENT_SCOPE)
