/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#include "gromacs/math/coordinatetransformation.h"

#include <array>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

using ::testing::Pointwise;
class TranslateAndScaleTest : public ::testing::Test
{
protected:
    RVec              identityScale_       = { 1, 1, 1 };
    RVec              identityTranslation_ = { 0, 0, 0 };
    std::vector<RVec> testVectors_         = { { 0, 0, 0 },
                                       { 1, 0, 0 },
                                       { 0, -1, -1 },
                                       { 1e10, 1e1, 1e-2 },
                                       { 3, -6, 2.5 } };
};

TEST_F(TranslateAndScaleTest, identityTransformation)
{
    TranslateAndScale translateAndScale(identityScale_, identityTranslation_);
    auto              toBeTransformed = testVectors_;
    translateAndScale(toBeTransformed);
    EXPECT_THAT(testVectors_, Pointwise(RVecEq(defaultFloatTolerance()), toBeTransformed));
}

TEST_F(TranslateAndScaleTest, translationWithIdentityScaling)
{
    TranslateAndScale translateAndScale(identityScale_, { 1, -1, 2 });
    translateAndScale(testVectors_);
    std::vector<RVec> expected = {
        { 1, -1, 2 }, { 2, -1, 2 }, { 1, -2, 1 }, { 1e10 + 1, 1e1 - 1, 1e-2 + 2 }, { 4, -7, 4.5 }
    };
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

TEST_F(TranslateAndScaleTest, scalingWithZeroTranslation)
{
    TranslateAndScale translateAndScale({ -1e10, 2, 0 }, identityTranslation_);
    translateAndScale(testVectors_);
    std::vector<RVec> expected = {
        { 0, 0, 0 }, { -1e10, 0, 0 }, { 0, -2, 0 }, { -1e20, 2e1, 0 }, { -3e10, -12, 0 }
    };
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

TEST_F(TranslateAndScaleTest, translationAndScalingNonTrivial)
{
    TranslateAndScale translateAndScale({ -1e10, 2, 0 }, { 1, -1, 2 });
    translateAndScale(testVectors_);
    std::vector<RVec> expected = {
        { -1e+10, -2, 0 }, { -2e+10, -2, 0 }, { -1e+10, -4, 0 }, { -1e+20, 18, 0 }, { -4e+10, -14, 0 }
    };
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

TEST_F(TranslateAndScaleTest, scalingIdentity)
{
    ScaleCoordinates scale(identityScale_);
    auto             expected = testVectors_;
    scale(testVectors_);
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

TEST_F(TranslateAndScaleTest, scalingNonTrivial)
{
    ScaleCoordinates  scale({ -100, 0.1, 0 });
    std::vector<RVec> expected = {
        { 0, 0, 0 }, { -100, 0, 0 }, { 0, -0.1, 0 }, { -1e12, 1e0, 0 }, { -300, -0.6, 0 }
    };
    scale(testVectors_);
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

TEST_F(TranslateAndScaleTest, scalingInverseNoZero)
{
    ScaleCoordinates scale({ -100, 0.1, 3 });
    auto             expected = testVectors_;
    scale(testVectors_);
    scale.inverseIgnoringZeroScale(testVectors_);
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

TEST_F(TranslateAndScaleTest, scalingInverseWithOneScaleDimensionZero)
{
    ScaleCoordinates scale({ -100, 0.1, 0 });
    std::vector<RVec> expected = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 1e10, 1e1, 0 }, { 3, -6, 0 } };
    scale(testVectors_);
    scale.inverseIgnoringZeroScale(testVectors_);
    EXPECT_THAT(expected, Pointwise(RVecEq(defaultFloatTolerance()), testVectors_));
}

} // namespace test
} // namespace gmx
