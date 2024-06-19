/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Tests \Gromacs-specific test assertions.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/testasserts.h"

#include <limits>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

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
template<typename FloatType>
FloatType addUlps(FloatType value, int ulps)
{
    return FloatingPoint<FloatType>::ReinterpretBits(FloatingPoint<FloatType>(value).bits() + ulps);
}

using ::gmx::test::FloatingPointDifference;

TEST(FloatingPointDifferenceTest, HandlesEqualValues)
{
    FloatingPointDifference diff(1.2, 1.2);
    EXPECT_TRUE(diff.isDouble());
    EXPECT_FALSE(diff.isNaN());
    EXPECT_EQ(0.0, diff.asAbsolute());
    EXPECT_EQ(0U, diff.asUlps());
    EXPECT_FALSE(diff.signsDiffer());
}

// TODO: Use typed tests to run all the tests for single and double.
TEST(FloatingPointDifferenceTest, HandlesFloatValues)
{
    FloatingPointDifference diff(1.2F, 1.2F);
    EXPECT_FALSE(diff.isDouble());
    EXPECT_FALSE(diff.isNaN());
    EXPECT_EQ(0.0, diff.asAbsolute());
    EXPECT_EQ(0U, diff.asUlps());
    EXPECT_FALSE(diff.signsDiffer());
}

TEST(FloatingPointDifferenceTest, HandlesZerosOfDifferentSign)
{
    FloatingPointDifference diff(0.0, GMX_DOUBLE_NEGZERO);
    EXPECT_FALSE(diff.isNaN());
    EXPECT_EQ(0.0, diff.asAbsolute());
    EXPECT_EQ(0U, diff.asUlps());
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
    const double            second = addUlps(0.0, 2);
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

TEST(FloatingPointToleranceTest, UlpTolerance)
{
    using gmx::test::ulpTolerance;

    FloatingPointDifference fequal(1.0, 1.0);
    FloatingPointDifference fulp2(1.0F, addUlps(1.0F, 2));
    EXPECT_TRUE(ulpTolerance(0).isWithin(fequal));
    EXPECT_FALSE(ulpTolerance(1).isWithin(fulp2));
    EXPECT_TRUE(ulpTolerance(2).isWithin(fulp2));

    FloatingPointDifference dequal(1.0, 1.0);
    FloatingPointDifference dulp2(1.0, addUlps(1.0, 2));
    FloatingPointDifference dulp2f(1.0, static_cast<double>(addUlps(1.0F, 2)));
    EXPECT_TRUE(ulpTolerance(0).isWithin(dequal));
    EXPECT_TRUE(ulpTolerance(2).isWithin(dulp2));
    EXPECT_FALSE(ulpTolerance(2).isWithin(dulp2f));
}

TEST(FloatingPointToleranceTest, RelativeToleranceAsFloatingPoint)
{
    using gmx::test::relativeToleranceAsFloatingPoint;

    FloatingPointDifference fequal(1.0F, 1.0F);
    FloatingPointDifference fulp2(1.0F, addUlps(1.0F, 2));
    FloatingPointDifference fdiff(1.0F, 1.011F);
    FloatingPointDifference fsmall(0.1F, 0.111F);
    FloatingPointDifference fsmall2(0.1F, 0.121F);
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 1e-2).isWithin(fequal));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 1e-9).isWithin(fequal));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 1e-2).isWithin(fulp2));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(1.0, 1e-9).isWithin(fulp2));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(1.0, 1e-2).isWithin(fdiff));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 2e-2).isWithin(fdiff));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 2e-2).isWithin(fsmall));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(0.1, 2e-2).isWithin(fsmall));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(1.0, 2e-2).isWithin(fsmall2));

    FloatingPointDifference dequal(1.0, 1.0);
    FloatingPointDifference dulp2f(1.0, static_cast<double>(addUlps(1.0F, 2)));
    FloatingPointDifference ddiff(1.0, 1.011);
    FloatingPointDifference dsmall(0.1, 0.111);
    FloatingPointDifference dsmall2(0.1, 0.121);
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 1e-2).isWithin(dequal));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 1e-9).isWithin(dequal));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 1e-2).isWithin(dulp2f));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(1.0, 1e-9).isWithin(dulp2f));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(1.0, 1e-2).isWithin(ddiff));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 2e-2).isWithin(ddiff));
    EXPECT_TRUE(relativeToleranceAsFloatingPoint(1.0, 2e-2).isWithin(dsmall));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(0.1, 2e-2).isWithin(dsmall));
    EXPECT_FALSE(relativeToleranceAsFloatingPoint(1.0, 2e-2).isWithin(dsmall2));
}

TEST(FloatingPointToleranceTest, RelativeToleranceAsUlp)
{
    using gmx::test::relativeToleranceAsUlp;

    FloatingPointDifference fequal(1.0F, 1.0F);
    FloatingPointDifference fulp4(1.0F, addUlps(1.0F, 4));
    FloatingPointDifference fsmall(0.1F, addUlps(1.0F, 2) - 0.9F);
    FloatingPointDifference fsmall2(0.1F, addUlps(1.0F, 6) - 0.9F);
    EXPECT_TRUE(relativeToleranceAsUlp(1.0, 2).isWithin(fequal));
    EXPECT_FALSE(relativeToleranceAsUlp(1.0, 2).isWithin(fulp4));
    EXPECT_TRUE(relativeToleranceAsUlp(1.0, 4).isWithin(fulp4));
    EXPECT_TRUE(relativeToleranceAsUlp(1.0, 4).isWithin(fsmall));
    EXPECT_FALSE(relativeToleranceAsUlp(0.1, 4).isWithin(fsmall));
    EXPECT_FALSE(relativeToleranceAsUlp(1.0, 4).isWithin(fsmall2));

    FloatingPointDifference dequal(1.0, 1.0);
    FloatingPointDifference dulp4(1.0, addUlps(1.0, 4));
    FloatingPointDifference dulp4f(1.0, static_cast<double>(addUlps(1.0F, 4)));
    FloatingPointDifference dsmall(0.1, addUlps(1.0, 2) - 0.9);
    FloatingPointDifference dsmall2(0.1, addUlps(1.0, 6) - 0.9);
    EXPECT_TRUE(relativeToleranceAsUlp(1.0, 2).isWithin(dequal));
    EXPECT_FALSE(relativeToleranceAsUlp(1.0, 2).isWithin(dulp4));
    EXPECT_TRUE(relativeToleranceAsUlp(1.0, 4).isWithin(dulp4));
    EXPECT_FALSE(relativeToleranceAsUlp(1.0, 4).isWithin(dulp4f));
    EXPECT_TRUE(relativeToleranceAsUlp(1.0, 4).isWithin(dsmall));
    EXPECT_FALSE(relativeToleranceAsUlp(0.1, 4).isWithin(dsmall));
    EXPECT_FALSE(relativeToleranceAsUlp(1.0, 4).isWithin(dsmall2));
}

TEST(FloatingPointToleranceTest, DefaultFloatTolerance)
{
    using gmx::test::defaultFloatTolerance;

    // Differences within 4 single-precision ULPs are within the tolerance
    FloatingPointDifference fequal(1.0F, 1.0F);
    FloatingPointDifference fulp4(1.0F, addUlps(1.0F, 4));
    FloatingPointDifference fulp8(1.0F, addUlps(1.0F, 8));
    FloatingPointDifference fsmall(0.1F, addUlps(1.0F, 2) - 0.9F);
    FloatingPointDifference fsmall2(0.1F, addUlps(1.0F, 6) - 0.9F);
    EXPECT_TRUE(defaultFloatTolerance().isWithin(fequal));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(fulp4));
    EXPECT_FALSE(defaultFloatTolerance().isWithin(fulp8));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(fsmall));
    EXPECT_FALSE(defaultFloatTolerance().isWithin(fsmall2));

    // Differences within 4 single-precision ULPs are still within the
    // tolerance, even when expressed as double-precision values.
    FloatingPointDifference dequal(1.0, 1.0);
    FloatingPointDifference dulp4(1.0, addUlps(1.0, 4));
    FloatingPointDifference dulp8(1.0, addUlps(1.0, 8));
    FloatingPointDifference dulp4f(1.0, static_cast<double>(addUlps(1.0F, 4)));
    FloatingPointDifference dulp8f(1.0, static_cast<double>(addUlps(1.0F, 8)));
    FloatingPointDifference dsmallf(0.1, static_cast<double>(addUlps(1.0F, 2) - 0.9F));
    FloatingPointDifference dsmall2f(0.1, static_cast<double>(addUlps(1.0F, 6) - 0.9F));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(dequal));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(dulp4));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(dulp8));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(dulp4f));
    EXPECT_FALSE(defaultFloatTolerance().isWithin(dulp8f));
    EXPECT_TRUE(defaultFloatTolerance().isWithin(dsmallf));
    EXPECT_FALSE(defaultFloatTolerance().isWithin(dsmall2f));
}

} // namespace
