/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests for simple math functions.
 *
 * Theoretically we'd like all math functions to match exactly, even in
 * floating-point, but that is not realistic without excessive corrections to
 * get the last ulp right - so in practice we target 1-2 ulp, which is VERY
 * tight. However, to avoid worrying users we are a bit more generous and
 * use 5-10 ulp as tolerance for the actual tests users will run through.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/functions.h"

#include <cmath>
#include <cstddef>
#include <cstdint>

#include <initializer_list>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(FunctionTest, StaticLog2)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<int>                result(11);

    // This needs to be expanded manually since it is evaluated at compile time,
    // and the compiler chokes if we put it as formal arguments to push_back
    result[0]  = gmx::StaticLog2<1>::value;
    result[1]  = gmx::StaticLog2<2>::value;
    result[2]  = gmx::StaticLog2<3>::value;
    result[3]  = gmx::StaticLog2<4>::value;
    result[4]  = gmx::StaticLog2<5>::value;
    result[5]  = gmx::StaticLog2<6>::value;
    result[6]  = gmx::StaticLog2<7>::value;
    result[7]  = gmx::StaticLog2<8>::value;
    result[8]  = gmx::StaticLog2<0xFFFFFFFF>::value;
    result[9]  = gmx::StaticLog2<9876543210>::value; // > 32 bits
    result[10] = gmx::StaticLog2<0xFFFFFFFFFFFFFFFFULL>::value;

    checker.checkSequence(result.begin(), result.end(), "StaticLog2");
}

TEST(FunctionTest, Log2I32Bit)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<int>                result;

    for (std::uint32_t i = 1; i <= 0xF; i++)
    {
        result.push_back(gmx::log2I(i));
    }

    for (std::uint32_t i = 0; i <= 0xF; i++)
    {
        result.push_back(gmx::log2I(static_cast<std::uint32_t>(0xFFFFFFF0 + i)));
    }
    checker.checkSequence(result.begin(), result.end(), "Log2I32Bit");
}

TEST(FunctionTest, Log2I64Bit)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<int>                result;

    for (std::uint64_t i = 1; i <= 0xF; i++)
    {
        result.push_back(gmx::log2I(i));
    }

    for (std::uint64_t i = 0; i <= 0x1F; i++)
    {
        result.push_back(gmx::log2I(static_cast<std::uint64_t>(0xFFFFFFF0ULL + i)));
    }

    for (std::uint64_t i = 0; i <= 0xF; i++)
    {
        result.push_back(gmx::log2I(static_cast<std::uint64_t>(0xFFFFFFFFFFFFFFF0ULL + i)));
    }
    checker.checkSequence(result.begin(), result.end(), "Log2I64Bit");
}

TEST(FunctionTest, GreatestCommonDivisor)
{
    EXPECT_EQ(8, gmx::greatestCommonDivisor(24, 64));
    EXPECT_EQ(8, gmx::greatestCommonDivisor(64, 24));
    EXPECT_EQ(1, gmx::greatestCommonDivisor(169, 289));
}

TEST(FunctionTest, InvsqrtFloat)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<float>              result;

    for (float f = 1.0; f < 10.0; f += 1.0)
    {
        result.push_back(gmx::invsqrt(f));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvsqrtFloat");
}

TEST(FunctionTest, InvsqrtDouble)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (double f = 1.0; f < 10.0; f += 1.0) // NOLINT(clang-analyzer-security.FloatLoopCounter)
    {
        result.push_back(gmx::invsqrt(f));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvsqrtDouble");
}

TEST(FunctionTest, InvsqrtInteger)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (int i = 1; i < 10; i++)
    {
        result.push_back(gmx::invsqrt(i));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvsqrtInteger");
}

TEST(FunctionTest, InvcbrtFloat)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<float>              result;

    for (int f : { -5, -4, -3, -2, -1, 1, 2, 3, 4 })
    {
        result.push_back(gmx::invcbrt(f));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvcbrtFloat");
}

TEST(FunctionTest, InvcbrtDouble)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (double d : { -5, -4, -3, -2, -1, 1, 2, 3, 4 })
    {
        result.push_back(gmx::invcbrt(d));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvcbrtDouble");
}

TEST(FunctionTest, InvcbrtInteger)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (int i : { -5, -4, -3, -2, -1, 1, 2, 3, 4 })
    {
        result.push_back(gmx::invcbrt(i));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvcbrtInteger");
}


TEST(FunctionTest, SixthrootFloat)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<float>              result;

    for (float f = 0; f < 10.0; f += 1.0)
    {
        result.push_back(gmx::sixthroot(f));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "SixthrootFloat");
}

TEST(FunctionTest, SixthrootDouble)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (double d = 0; d < 10.0; d += 1.0) // NOLINT(clang-analyzer-security.FloatLoopCounter)
    {
        result.push_back(gmx::sixthroot(d));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "SixthrootDouble");
}

TEST(FunctionTest, SixthrootInteger)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    result.reserve(10);
    for (int i = 0; i < 10; i++)
    {
        result.push_back(gmx::sixthroot(i));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "SixthrootInteger");
}

TEST(FunctionTest, InvsixthrootFloat)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<float>              result;

    for (float f = 1.0; f < 10.0; f += 1.0)
    {
        result.push_back(gmx::invsixthroot(f));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvsixthrootFloat");
}

TEST(FunctionTest, InvsixthrootDouble)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (double d = 1.0; d < 10.0; d += 1.0) // NOLINT(clang-analyzer-security.FloatLoopCounter)
    {
        result.push_back(gmx::invsixthroot(d));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvsixthrootDouble");
}

TEST(FunctionTest, InvsixthrootInteger)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;

    for (int i = 1; i < 10; i++)
    {
        result.push_back(gmx::invsixthroot(i));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "InvsixthrootInteger");
}

TEST(FunctionTest, Powers)
{
    // These should be remarkably difficult to screw up, but test each
    // of them once anyway with integer and fp arguments to catch typos.
    EXPECT_EQ(4, gmx::square(2));
    EXPECT_EQ(8, gmx::power3(2));
    EXPECT_EQ(16, gmx::power4(2));
    EXPECT_EQ(32, gmx::power5(2));
    EXPECT_EQ(64, gmx::power6(2));
    EXPECT_EQ(4096, gmx::power12(2));

    EXPECT_EQ(6.25, gmx::square(2.5));
    EXPECT_EQ(15.625, gmx::power3(2.5));
    EXPECT_EQ(39.0625, gmx::power4(2.5));
    EXPECT_EQ(97.65625, gmx::power5(2.5));
    EXPECT_EQ(244.140625, gmx::power6(2.5));
    EXPECT_EQ(59604.644775390625, gmx::power12(2.5));
}

TEST(FunctionTest, ErfInvFloat)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<float>              result;
    int                             npoints = 10;

    for (int i = 0; i < npoints; i++)
    {
        float r = float(2 * i - npoints + 1) / float(npoints);

        result.push_back(gmx::erfinv(r));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "ErfInvFloat");
}

TEST(FunctionTest, ErfInvDouble)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    std::vector<double>             result;
    int                             npoints = 10;

    for (int i = 0; i < npoints; i++)
    {
        double r = double(2 * i - npoints + 1) / npoints;
        result.push_back(gmx::erfinv(r));
    }
    checker.setDefaultTolerance(gmx::test::ulpTolerance(5));
    checker.checkSequence(result.begin(), result.end(), "ErfInvDouble");
}

TEST(FunctionTest, ErfAndErfInvAreInversesFloat)
{
    int npoints = 1000;

    for (int i = 0; i < npoints; i++)
    {
        float r = float(2 * i - npoints + 1) / float(npoints);
        // There is extra rounding when we do two function calls,
        // and since the derivative is >1 we lose a few ulp extra
        // when going back to r.
        EXPECT_FLOAT_EQ_TOL(r, std::erf(gmx::erfinv(r)), gmx::test::ulpTolerance(10));
    }
}

TEST(FunctionTest, ErfAndErfInvAreInversesDouble)
{
    int npoints = 1000;

    for (int i = 0; i < npoints; i++)
    {
        double r = double(2 * i - npoints + 1) / npoints;
        // There is extra rounding when we do two function calls,
        // and since the derivative is >1 we lose a few ulp extra
        // when going back to r.
        EXPECT_DOUBLE_EQ_TOL(r, std::erf(gmx::erfinv(r)), gmx::test::ulpTolerance(10));
    }
}

template<typename T>
class FunctionTestIntegerTypes : public ::testing::Test
{
};

typedef ::testing::Types<signed char, unsigned char, short, unsigned short, int, unsigned int, long, unsigned long> IntegerTypes;
TYPED_TEST_SUITE(FunctionTestIntegerTypes, IntegerTypes);

TYPED_TEST(FunctionTestIntegerTypes, IsPowerOfTwo)
{
    if (std::is_signed_v<TypeParam>)
    {
        EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(std::numeric_limits<TypeParam>::min()));
        EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(-16));
        EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(-3));
        EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(-2));
        EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(-1));
    }
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(0));
    EXPECT_EQ(true, gmx::isPowerOfTwo<TypeParam>(1));
    EXPECT_EQ(true, gmx::isPowerOfTwo<TypeParam>(2));
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(3));
    EXPECT_EQ(true, gmx::isPowerOfTwo<TypeParam>(4));
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(5));
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(6));
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(24));
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(63));
    EXPECT_EQ(true, gmx::isPowerOfTwo<TypeParam>(64));
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(66));
    // Max for any type is always 2^x - 1
    EXPECT_EQ(false, gmx::isPowerOfTwo<TypeParam>(std::numeric_limits<TypeParam>::max()));
}


TYPED_TEST(FunctionTestIntegerTypes, DivideRoundUp)
{
    // Test some reasonable values and edge-cases
    // Constraint: the sum of numerator and denominator must fit into TypeParam
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(0, 1), 0);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(0, 10), 0);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(0, std::numeric_limits<TypeParam>::max()), 0);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(1, 1), 1);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(10, 1), 10);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(std::numeric_limits<TypeParam>::max() - 1, 1),
              std::numeric_limits<TypeParam>::max() - 1);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(2, 2), 1);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(3, 2), 2);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(4, 2), 2);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(5, 2), 3);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(1, 10), 1);
    EXPECT_EQ(gmx::divideRoundUp<TypeParam>(1, std::numeric_limits<TypeParam>::max() - 1), 1);

    // Test random inputs; up to square root of max value
    auto                                       rng = std::minstd_rand{ 20240207 };
    std::uniform_int_distribution<std::size_t> distrib(
            0, std::sqrt(std::numeric_limits<TypeParam>::max()));
    for (int iter = 0; iter < 400; iter++)
    {
        TypeParam nom    = distrib(rng);
        TypeParam denom  = distrib(rng) + 1;
        TypeParam result = gmx::divideRoundUp(nom, denom);
        EXPECT_GT(nom, denom * (result - 1));
        EXPECT_LE(nom, denom * result);
    }
}

} // namespace
} // namespace test
} // namespace gmx
