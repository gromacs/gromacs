/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Tests \Gromacs-specific test assertions.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "testutils/testasserts.h"

#include <gtest/gtest.h>

namespace
{

using ::testing::internal::FloatingPoint;

/*! \brief
 * Helper to produce floating-point numbers with specified ULP difference.
 *
 * This doesn't work if the value would change sign.
 *
 * \ingroup module_testutils
 */
template <typename FloatType>
FloatType addUlps(FloatType value, int ulps)
{
    return FloatingPoint<FloatType>::ReinterpretBits(
            FloatingPoint<FloatType>(value).bits() + ulps);
}

using ::gmx::test::FloatingPointDifference;

TEST(FloatingPointDifferenceTest, HandlesEqualValues)
{
    FloatingPointDifference diff(1.2, 1.2);
    EXPECT_FALSE(diff.isNaN());
    EXPECT_EQ(0.0, diff.asAbsolute());
    EXPECT_EQ(0U,  diff.asUlps());
    EXPECT_FALSE(diff.signsDiffer());
}

TEST(FloatingPointDifferenceTest, HandlesZerosOfDifferentSign)
{
    FloatingPointDifference diff(0.0, GMX_DOUBLE_NEGZERO);
    EXPECT_FALSE(diff.isNaN());
    EXPECT_EQ(0.0, diff.asAbsolute());
    EXPECT_EQ(0U,  diff.asUlps());
    EXPECT_TRUE(diff.signsDiffer());
}

TEST(FloatingPointDifferenceTest, HandlesSignComparisonWithZero)
{
    {
        FloatingPointDifference diff(0.0, -1.2);
        EXPECT_FALSE(diff.isNaN());
        EXPECT_DOUBLE_EQ(1.2, diff.asAbsolute());
        EXPECT_TRUE(diff.signsDiffer());
    }
    {
        FloatingPointDifference diff(GMX_DOUBLE_NEGZERO, -1.2);
        EXPECT_FALSE(diff.isNaN());
        EXPECT_DOUBLE_EQ(1.2, diff.asAbsolute());
        EXPECT_FALSE(diff.signsDiffer());
    }
}

TEST(FloatingPointDifferenceTest, HandlesUlpDifferences)
{
    const double first  = 1.0;
    const double second = addUlps(first, 2);
    {
        FloatingPointDifference diff(first, second);
        EXPECT_FALSE(diff.isNaN());
        EXPECT_DOUBLE_EQ(second - first, diff.asAbsolute());
        EXPECT_EQ(2U, diff.asUlps());
        EXPECT_FALSE(diff.signsDiffer());
    }
    {
        FloatingPointDifference diff(second, first);
        EXPECT_FALSE(diff.isNaN());
        EXPECT_DOUBLE_EQ(second - first, diff.asAbsolute());
        EXPECT_EQ(2U, diff.asUlps());
        EXPECT_FALSE(diff.signsDiffer());
    }
}

TEST(FloatingPointDifferenceTest, HandlesUlpDifferenceAcrossZero)
{
    const double            first  = addUlps(GMX_DOUBLE_NEGZERO, 2);
    const double            second = addUlps( 0.0, 2);
    FloatingPointDifference diff(first, second);
    EXPECT_FALSE(diff.isNaN());
    EXPECT_DOUBLE_EQ(second - first, diff.asAbsolute());
    EXPECT_EQ(4U, diff.asUlps());
    EXPECT_TRUE(diff.signsDiffer());
}

TEST(FloatingPointDifferenceTest, HandlesNaN)
{
    const double            first  = std::numeric_limits<double>::quiet_NaN();
    const double            second = 2.0;
    FloatingPointDifference diff(first, second);
    EXPECT_TRUE(diff.isNaN());
}

} // namespace
