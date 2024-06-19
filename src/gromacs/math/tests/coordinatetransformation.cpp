/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests structure similarity measures rmsd and size-independent rho factor.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/coordinatetransformation.h"

#include <array>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/matrix.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

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

class AffineTransformationTest : public ::testing::Test
{
protected:
    std::vector<RVec> testVectors_ = { { 0, 0, 0 },
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

TEST_F(TranslateAndScaleTest, translationAndScalingNonTrivialSingeVector)
{
    TranslateAndScale translateAndScale({ -1e10, 2, 0 }, { 1, -1, 2 });
    RVec              test(0, 0, 0);
    translateAndScale(&test);

    EXPECT_REAL_EQ(-1e+10, test[XX]);
    EXPECT_REAL_EQ(-2, test[YY]);
    EXPECT_REAL_EQ(0, test[ZZ]);
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

TEST_F(TranslateAndScaleTest, scalingNonTrivialSingleVector)
{
    ScaleCoordinates scale({ -100, 0.1, 0 });
    RVec             test(3, -6, 2.5);
    scale(&test);

    EXPECT_REAL_EQ(-300, test[XX]);
    EXPECT_REAL_EQ(-0.6, test[YY]);
    EXPECT_REAL_EQ(0, test[ZZ]);
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

TEST_F(TranslateAndScaleTest, scalingInverseWithOneScaleDimensionZeroSingleVector)
{
    ScaleCoordinates scale({ -100, 0.1, 0 });
    RVec             test(3, -6, 2.5);

    scale(&test);
    scale.inverseIgnoringZeroScale(&test);

    EXPECT_REAL_EQ(3, test[XX]);
    EXPECT_REAL_EQ(-6, test[YY]);
    EXPECT_REAL_EQ(0, test[ZZ]);
}

TEST_F(AffineTransformationTest, identityTransformYieldsSameVectors)
{
    const AffineTransformation identityTransformation(identityMatrix<real, 3>(), { 0, 0, 0 });
    for (const auto& vector : testVectors_)
    {
        RVec vectorTransformed = vector;
        identityTransformation(&vectorTransformed);
        EXPECT_REAL_EQ(vector[XX], vectorTransformed[XX]);
        EXPECT_REAL_EQ(vector[YY], vectorTransformed[YY]);
        EXPECT_REAL_EQ(vector[ZZ], vectorTransformed[ZZ]);
    }
}

TEST_F(AffineTransformationTest, applyTransformationToVectors)
{
    const Matrix3x3 transformMatrix({ 0.1, 1, 0.1, 0.4, 1, 0.6, 0.7, 0.8, 0.9 });
    const RVec      transformVector = { 1, -1e5, 1e4 };

    const AffineTransformation affineTransformation(transformMatrix, transformVector);

    const std::vector<RVec> expectedResult = { { 1, -100'000, 10'000 },
                                               { 1.1, -99999.6, 10000.7 },
                                               { -0.1, -100'002, 9998.3 },
                                               { 1e9, 3.9999e9, 7.00001e9 },
                                               { -4.45, -100'003, 9'999.5 } };

    auto expected = expectedResult.begin();
    for (const auto& vector : testVectors_)
    {
        RVec vectorTransformed = vector;
        affineTransformation(&vectorTransformed);
        // need relaxed tolerance here, due to the number of operations involved
        EXPECT_REAL_EQ_TOL((*expected)[XX],
                           vectorTransformed[XX],
                           relativeToleranceAsFloatingPoint((*expected)[XX], 1e-5));
        EXPECT_REAL_EQ_TOL((*expected)[YY],
                           vectorTransformed[YY],
                           relativeToleranceAsFloatingPoint((*expected)[YY], 1e-5));
        EXPECT_REAL_EQ_TOL((*expected)[ZZ],
                           vectorTransformed[ZZ],
                           relativeToleranceAsFloatingPoint((*expected)[ZZ], 1e-5));
        ++expected;
    }
}

TEST_F(AffineTransformationTest, retrieveGradient)
{
    const Matrix3x3            transformMatrix({ 0.1, 1, 0.1, 0.4, 1, 0.6, 0.7, 0.8, 0.9 });
    const RVec                 transformVector = { 1, -1e5, 1e4 };
    const AffineTransformation affineTransformation(transformMatrix, transformVector);

    const Matrix3x3 gradient = affineTransformation.gradient();

    const Matrix3x3 expectedResult({ 0.1, 0.4, 0.7, 1, 1, 0.8, 0.1, 0.6, 0.9 });

    for (int row = 0; row < 3; row++)
    {
        for (int column = 0; column < 3; column++)
        {
            EXPECT_REAL_EQ(gradient(row, column), expectedResult(row, column));
        }
    }
}

} // namespace test
} // namespace gmx
