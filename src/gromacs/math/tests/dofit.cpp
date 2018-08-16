/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Tests structure similarity measures rmsd and size-independent rho factor.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include <array>

#include <gtest/gtest.h>

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace
{

using ::testing::Pointwise;
using gmx::RVec;
using gmx::test::defaultRealTolerance;
using gmx::findStructureSimilarity;
using gmx::RMSD;
using gmx::RhoMeasure;
class StructureSimilarityTest : public ::testing::Test
{
    protected:
        static constexpr int     nAtoms = 4;
        std::array<RVec, nAtoms> structureA {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};
        std::array<RVec, nAtoms> structureB {{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {0, 0, 0}}};
        std::array<real, nAtoms> masses     {{1, 1, 1, 0}};
        std::array<int, 3>       index      {{0, 1, 2}};
        rvec                   * x1  = gmx::as_rvec_array(structureA.data());
        rvec                   * x2  = gmx::as_rvec_array(structureB.data());
        real                   * m   = masses.data();
};

TEST_F(StructureSimilarityTest, StructureComparedToSelfHasZeroRMSD)
{
    EXPECT_REAL_EQ_TOL(0., findStructureSimilarity<RMSD>(nAtoms, m, x1, x1), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, StructureComparedToSelfHasZeroRho)
{
    EXPECT_REAL_EQ_TOL(0., findStructureSimilarity<RhoMeasure>(nAtoms, m, x1, x1), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRMSD)
{
    EXPECT_REAL_EQ_TOL(sqrt(2.0), findStructureSimilarity<RMSD>(nAtoms, m, x1, x2), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRho)
{
    EXPECT_REAL_EQ_TOL(2., findStructureSimilarity<RhoMeasure>(nAtoms, m, x1, x2), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRMSDWithIndex)
{
    EXPECT_REAL_EQ_TOL(sqrt(2.0), findStructureSimilarity<RMSD>(index.size(), m, x1, x2, index.data()), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRhoWidthIndex)
{
    EXPECT_REAL_EQ_TOL(2., findStructureSimilarity<RhoMeasure>(index.size(), m, x1, x2, index.data()), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectFitRotationMatrixThreeD)
{
    constexpr matrix expected = {{0, 1, 0},
                                 {0, 0, 1},
                                 {1, 0, 0}};
    matrix           R;
    calc_fit_R(3, nAtoms, m, x1, x2, R);
    EXPECT_THAT(R, Pointwise(RVecEq(defaultRealTolerance()), expected));
}

TEST_F(StructureSimilarityTest, YieldsCorrectFitRotationMatrixTwoD)
{
    constexpr matrix expected = {{0, 1, 0},
                                 {-1, 0, 0},
                                 {0, 0, 1}};
    matrix           R;
    calc_fit_R(2, nAtoms, m, x1, x2, R);
    EXPECT_THAT(R, Pointwise(RVecEq(defaultRealTolerance()), expected));
}

} // namespace
