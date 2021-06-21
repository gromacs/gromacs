/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <numeric>

#include <gtest/gtest.h>

#include "nblib/listed_forces/conversions.hpp"

#include "testutils/testasserts.h"


namespace nblib
{
namespace test
{
namespace
{

ListedInteractionData someBondsAndAngles()
{
    ListedInteractionData         interactions;
    HarmonicBondType              bond1{ 10, 0.1 };
    HarmonicBondType              bond2{ 20, 0.2 };
    std::vector<HarmonicBondType> bonds{ bond1, bond2 };
    pickType<HarmonicBondType>(interactions).parameters = bonds;

    HarmonicAngle              angle1(100, Degrees(100));
    HarmonicAngle              angle2(200, Degrees(101));
    std::vector<HarmonicAngle> angles{ angle1, angle2 };
    pickType<HarmonicAngle>(interactions).parameters = angles;

    std::vector<InteractionIndex<HarmonicBondType>> bondIndices{ { 0, 1, 0 }, { 1, 2, 0 }, { 2, 3, 1 } };
    pickType<HarmonicBondType>(interactions).indices = std::move(bondIndices);

    std::vector<InteractionIndex<HarmonicAngle>> angleIndices{ { 0, 1, 2, 0 }, { 1, 2, 3, 1 } };
    pickType<HarmonicAngle>(interactions).indices = std::move(angleIndices);

    return interactions;
}

TEST(ListedShims, ParameterConversion)
{
    ListedInteractionData interactions = someBondsAndAngles();

    auto [idef, gmx_params] = createFFparams(interactions);

    EXPECT_EQ(gmx_params->iparams.size(), 4);
    EXPECT_EQ(gmx_params->iparams[0].harmonic.rA,
              pickType<HarmonicBondType>(interactions).parameters[0].equilConstant());
    EXPECT_REAL_EQ_TOL(gmx_params->iparams[2].harmonic.rA,
                       pickType<HarmonicAngle>(interactions).parameters[0].equilConstant() / DEG2RAD,
                       gmx::test::defaultRealTolerance());

    EXPECT_EQ(idef->il[F_BONDS].iatoms.size(), 9);
    std::vector<int> bondIatoms{ 0, 0, 1, 0, 1, 2, 1, 2, 3 };
    EXPECT_EQ(idef->il[F_BONDS].iatoms, bondIatoms);
    std::vector<int> angleIatoms{ 2, 0, 1, 2, 3, 1, 2, 3 };
    EXPECT_EQ(idef->il[F_ANGLES].iatoms, angleIatoms);
    idef->clear();
}

} // namespace
} // namespace test
} // namespace nblib
