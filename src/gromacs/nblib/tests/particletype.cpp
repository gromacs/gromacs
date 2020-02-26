/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * This implements basic nblib AtomType tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 */
#include "gmxpre.h"

#include "gromacs/nblib/particletype.h"

#include <cmath>

#include "testutils/testasserts.h"

#include "testsystems.h"

using gmx::test::defaultRealTolerance;

namespace nblib
{

TEST(NBlibTest, ParticleTypeNameCanBeConstructed)
{
    ArAtom       arAtom;
    ParticleType argonAtom(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    EXPECT_EQ(argonAtom.name(), arAtom.name);
}

TEST(NBlibTest, ParticleTypeMassCanBeConstructed)
{
    ArAtom       arAtom;
    ParticleType argonAtom(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    EXPECT_EQ(argonAtom.mass(), arAtom.mass);
}

TEST(NBlibTest, ParticleTypeC6CanBeConstructed)
{
    ArAtom       arAtom;
    ParticleType argonAtom(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    EXPECT_EQ(argonAtom.c6(), arAtom.c6);
}

TEST(NBlibTest, ParticleTypeC12CanBeConstructed)
{
    ArAtom       arAtom;
    ParticleType argonAtom(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);
    EXPECT_EQ(argonAtom.c12(), arAtom.c12);
}

} // namespace nblib
