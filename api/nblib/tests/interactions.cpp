/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * This implements molecule setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/interactions.h"

#include <cmath>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

#include "nblib/exception.h"
#include "nblib/particletype.h"
#include "nblib/tests/testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, NonBondedForceParamsCorrect)
{
    ParticleType atom1(ParticleTypeName("a1"), Mass(1));
    ParticleType atom2(ParticleTypeName("a2"), Mass(1));
    ParticleType atom3(ParticleTypeName("a3"), Mass(1));

    ParticleTypesInteractions interactions;

    auto c6_1  = C6(1.6);
    auto c12_1 = C12(1.12);
    C6   c6_2{ 2.6 };
    C12  c12_2{ 2.12 };
    C6   c6_3  = C6(3.6);
    C12  c12_3 = C12(3.12);

    auto c6comb  = C6{ 40 };
    auto c12comb = C12{ 41 };

    interactions.add(atom1.name(), c6_1, c12_1);
    interactions.add(atom2.name(), c6_2, c12_2);
    interactions.add(atom3.name(), c6_3, c12_3);
    interactions.add(atom2.name(), atom3.name(), c6comb, c12comb);

    auto nbfp = interactions.generateTable();

    //! self interaction for c6
    EXPECT_REAL_EQ_TOL(c6_1, nbfp.getC6(atom1.name(), atom1.name()), gmx::test::defaultRealTolerance());
    //! + symmetric pair
    EXPECT_REAL_EQ_TOL(c12_1, nbfp.getC12(atom1.name(), atom1.name()), gmx::test::defaultRealTolerance());

    //! geometric comb rule for c6
    EXPECT_REAL_EQ_TOL(std::sqrt(c6_1 * c6_2),
                       nbfp.getC6(atom1.name(), atom2.name()),
                       gmx::test::defaultRealTolerance());
    //! + symmetric pair
    EXPECT_REAL_EQ_TOL(std::sqrt(c6_1 * c6_2),
                       nbfp.getC6(atom2.name(), atom1.name()),
                       gmx::test::defaultRealTolerance());

    //! geometric comb rule for c12
    EXPECT_REAL_EQ_TOL(std::sqrt(c12_1 * c12_2),
                       nbfp.getC12(atom1.name(), atom2.name()),
                       gmx::test::defaultRealTolerance());

    //! + symmetric par
    EXPECT_REAL_EQ_TOL(std::sqrt(c12_1 * c12_2),
                       nbfp.getC12(atom2.name(), atom1.name()),
                       gmx::test::defaultRealTolerance());

    //! explicit pairwise interaction c6
    EXPECT_REAL_EQ_TOL(c6comb, nbfp.getC6(atom2.name(), atom3.name()), gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(c6comb, nbfp.getC6(atom3.name(), atom2.name()), gmx::test::defaultRealTolerance());

    //! explicit pairwise interaction c12
    EXPECT_REAL_EQ_TOL(
            c12comb, nbfp.getC12(atom2.name(), atom3.name()), gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(
            c12comb, nbfp.getC12(atom3.name(), atom2.name()), gmx::test::defaultRealTolerance());

    ParticleType atom4(ParticleTypeName("a4"), Mass(1));
    interactions.add(atom3.name(), atom4.name(), C6(1), C12(2));
    //! need to have self-interaction for all particles
    EXPECT_THROW(auto interactionsTable = interactions.generateTable(), InputException);
}

TEST(NBlibTest, CanMergeInteractions)
{
    ParticleType atom1(ParticleTypeName("a1"), Mass(1));
    ParticleType atom2(ParticleTypeName("a2"), Mass(1));
    ParticleType atom3(ParticleTypeName("a3"), Mass(1));

    ParticleTypesInteractions interactions;

    auto c6_1  = C6(1.6);
    auto c12_1 = C12(1.12);
    C6   c6_2{ 2.6 };
    C12  c12_2{ 2.12 };
    C6   c6_3  = C6(3.6);
    C12  c12_3 = C12(3.12);

    auto c6comb  = C6{ 40 };
    auto c12comb = C12{ 41 };

    interactions.add(atom1.name(), c6_1, c12_1);
    interactions.add(atom2.name(), c6_2, c12_2);
    interactions.add(atom3.name(), c6_3, c12_3);
    interactions.add(atom2.name(), atom3.name(), c6comb, c12comb);

    ParticleType atom4(ParticleTypeName("a4"), Mass(1));
    ParticleType atom5(ParticleTypeName("a5"), Mass(1));
    auto         c6_4  = C6{ 4.6 };
    auto         c12_4 = C12{ 4.12 };

    C6  c6_5{ 5.6 };
    C12 c12_5{ 5.12 };

    C6  c6_override  = C6(45.6);
    C12 c12_override = C12(45.12);

    ParticleTypesInteractions otherInteractions;
    otherInteractions.add(atom4.name(), c6_4, c12_4);
    otherInteractions.add(atom5.name(), c6_5, c12_5);
    otherInteractions.add(atom4.name(), atom5.name(), c6_override, c12_override);

    interactions.merge(otherInteractions);

    auto nbfp = interactions.generateTable();

    EXPECT_REAL_EQ_TOL(std::sqrt(c6_3 * c6_4),
                       nbfp.getC6(atom3.name(), atom4.name()),
                       gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(std::sqrt(c12_3 * c12_4),
                       nbfp.getC12(atom3.name(), atom4.name()),
                       gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(
            c6_override, nbfp.getC6(atom4.name(), atom5.name()), gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(
            c12_override, nbfp.getC12(atom4.name(), atom5.name()), gmx::test::defaultRealTolerance());
}

} // namespace
} // namespace test
} // namespace nblib
