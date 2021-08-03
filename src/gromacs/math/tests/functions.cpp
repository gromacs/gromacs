/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018,2019 by the GROMACS development team.
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * Tests for simple math functions.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/functions.h"

#include <cmath>
#include <cstdint>

#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

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
    checker.checkSequence(result.begin(), result.end(), "ErfInvDouble");
}

TEST(FunctionTest, ErfAndErfInvAreInversesFloat)
{
    int npoints = 1000;

    for (int i = 0; i < npoints; i++)
    {
        float r = float(2 * i - npoints + 1) / float(npoints);
        EXPECT_FLOAT_EQ_TOL(r, std::erf(gmx::erfinv(r)), gmx::test::ulpTolerance(10));
    }
}

TEST(FunctionTest, ErfAndErfInvAreInversesDouble)
{
    int npoints = 1000;

    for (int i = 0; i < npoints; i++)
    {
        double r = double(2 * i - npoints + 1) / npoints;
        EXPECT_DOUBLE_EQ_TOL(r, std::erf(gmx::erfinv(r)), gmx::test::ulpTolerance(10));
    }
}

template<typename T>
class FunctionTestIntegerTypes : public ::testing::Test
{
};

typedef ::testing::Types<char, unsigned char, int, unsigned int, long, unsigned long> IntegerTypes;
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

} // namespace
