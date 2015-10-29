/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

TEST(FunctionTest, StaticLog2)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<int>                 result(10);

    // This needs to be expanded manually since it is evaluated at compile time,
    // and the compiler chokes if we put it as formal arguments to push_back
    result[0] = gmx::StaticLog2<1>::value;
    result[1] = gmx::StaticLog2<2>::value;
    result[2] = gmx::StaticLog2<3>::value;
    result[3] = gmx::StaticLog2<4>::value;
    result[4] = gmx::StaticLog2<5>::value;
    result[5] = gmx::StaticLog2<6>::value;
    result[6] = gmx::StaticLog2<7>::value;
    result[7] = gmx::StaticLog2<8>::value;
    result[8] = gmx::StaticLog2<9>::value;
    result[9] = gmx::StaticLog2<9876543210>::value; // > 32 bits

    checker.checkSequence(result.begin(), result.end(), "StaticLog2");

}

TEST(FunctionTest, GreatestCommonDivisor)
{
    EXPECT_EQ(8, gmx::greatestCommonDivisor(24, 64));
    EXPECT_EQ(8, gmx::greatestCommonDivisor(64, 24));
    EXPECT_EQ(1, gmx::greatestCommonDivisor(169, 289));
}

TEST(FunctionTest, InvsqrtFloat)
{
    // Compare our software invsqrt to the system 1/sqrt
    for (float x = 1e-3; x < 1000.0; x += 1e-3)
    {
        // Target is 22 bits in single, and then we add 2 bits margin (rounding errors).
        // Since the single precision fraction is 23 bits, this gives the ulp
        // tolerance (1 << (2 + 23 - 22) = 8
        EXPECT_REAL_EQ_TOL(1.0/sqrt(x), gmx::invsqrt(x), gmx::test::ulpTolerance(8));
    }
}

TEST(FunctionTest, InvsqrtDouble)
{
    // Compare our software invsqrt to the system 1/sqrt
    for (double x = 1e-3; x < 1000.0; x += 1e-3)
    {
        // Target is 44 bits in double, and then we add 2 bits margin (rounding errors).
        // Since the single precision fraction is 52 bits, this gives the ulp
        // tolerance (1 << (2 + 52 - 44) = 1024
        EXPECT_REAL_EQ_TOL(1.0/sqrt(x), gmx::invsqrt(x), gmx::test::ulpTolerance(1024));
    }
}

TEST(FunctionTest, InvsqrtInteger)
{
    // Compare our software invsqrt to the system 1/sqrt
    for (int i = 1; i < 1000; i += 1)
    {
        // Target is 44 bits in double, and then we add 2 bits margin (rounding errors).
        // Since the single precision fraction is 52 bits, this gives the ulp
        // tolerance (1 << (2 + 52 - 44) = 1024
        EXPECT_REAL_EQ_TOL(1.0/sqrt(i), gmx::invsqrt(i), gmx::test::ulpTolerance(1024));
    }
}

TEST(FunctionTest, FastsqrtFloat)
{
    // Compare our software fastsqrt to the system sqrt
    for (float x = 0; x < 1000.0; x += 1e-3)
    {
        // Target is 22 bits in single, and then we add 2 bits margin (rounding errors).
        // Since the single precision fraction is 23 bits, this gives the ulp
        // tolerance (1 << (2 + 23 - 22) = 8
        EXPECT_REAL_EQ_TOL(sqrt(x), gmx::fastsqrt(x), gmx::test::ulpTolerance(8));
    }
}

TEST(FunctionTest, FastsqrtDouble)
{
    // Compare our software fastsqrt to the system sqrt
    for (double x = 0; x < 1000.0; x += 1e-3)
    {
        // Target is 44 bits in double, and then we add 2 bits margin (rounding errors).
        // Since the single precision fraction is 52 bits, this gives the ulp
        // tolerance (1 << (2 + 52 - 44) = 1024
        EXPECT_REAL_EQ_TOL(sqrt(x), gmx::fastsqrt(x), gmx::test::ulpTolerance(1024));
    }
}

TEST(FunctionTest, FastsqrtInteger)
{
    // Compare our software fastsqrt to the system sqrt
    for (int i = 0; i < 1000; i += 1)
    {
        // Target is 44 bits in double, and then we add 2 bits margin (rounding errors).
        // Since the single precision fraction is 52 bits, this gives the ulp
        // tolerance (1 << (2 + 52 - 44) = 1024
        EXPECT_REAL_EQ_TOL(sqrt(i), gmx::fastsqrt(i), gmx::test::ulpTolerance(1024));
    }
}

TEST(FunctionTest, InvcbrtFloat)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<float>               result;

    for (float x = -5; x < 0; x += 1.0)
    {
        result.push_back(gmx::invcbrt(x));
    }
    for (float x = 1; x < 5; x += 1.0)
    {
        result.push_back(gmx::invcbrt(x));
    }

    checker.checkSequence(result.begin(), result.end(), "InvcbrtFloat");
}

TEST(FunctionTest, InvcbrtDouble)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<double>              result;

    for (double x = -5; x < 0; x += 1.0)
    {
        result.push_back(gmx::invcbrt(x));
    }
    for (double x = 1; x < 5; x += 1.0)
    {
        result.push_back(gmx::invcbrt(x));
    }
    checker.checkSequence(result.begin(), result.end(), "InvcbrtDouble");
}

TEST(FunctionTest, InvcbrtInteger)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<double>              result;

    for (int x = -5; x < 0; x += 1.0)
    {
        result.push_back(gmx::invcbrt(x));
    }
    for (int x = 1; x < 5; x += 1.0)
    {
        result.push_back(gmx::invcbrt(x));
    }
    checker.checkSequence(result.begin(), result.end(), "InvcbrtInteger");
}


TEST(FunctionTest, SixthrootFloat)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<float>               result;

    for (float x = 0; x < 100.0; x += 1.0)
    {
        result.push_back(gmx::sixthroot(x));
    }
    checker.checkSequence(result.begin(), result.end(), "SixthrootFloat");
}

TEST(FunctionTest, SixthrootDouble)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<double>              result;

    for (double x = 0; x < 100.0; x += 1.0)
    {
        result.push_back(gmx::sixthroot(x));
    }
    checker.checkSequence(result.begin(), result.end(), "SixthrootDouble");
}

TEST(FunctionTest, SixthrootInteger)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<double>              result;

    for (int i = 0; i < 100; i++)
    {
        result.push_back(gmx::sixthroot(i));
    }
    checker.checkSequence(result.begin(), result.end(), "SixthrootInteger");
}

TEST(FunctionTest, InvsixthrootFloat)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<float>               result;

    for (float x = 1e-3; x < 100.0; x += 1.0)
    {
        result.push_back(gmx::sixthroot(x));
    }
    checker.checkSequence(result.begin(), result.end(), "InvsixthrootFloat");
}

TEST(FunctionTest, InvsixthrootDouble)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<double>              result;

    for (double x = 1e-3; x < 100.0; x += 1.0)
    {
        result.push_back(gmx::sixthroot(x));
    }
    checker.checkSequence(result.begin(), result.end(), "InvsixthrootDouble");
}

TEST(FunctionTest, InvsixthrootInteger)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());
    std::vector<double>              result;

    for (int i = 1; i < 100; i++)
    {
        result.push_back(gmx::sixthroot(i));
    }
    checker.checkSequence(result.begin(), result.end(), "InvsixthrootInteger");
}

TEST(FunctionTest, Powers)
{
    // These should be remarkably difficult to screw up, but test each
    // of them once anyway with integer and fp arguments to catch typos.
    EXPECT_EQ(4,   gmx::square(2));
    EXPECT_EQ(8,   gmx::power3(2));
    EXPECT_EQ(16,  gmx::power4(2));
    EXPECT_EQ(32,  gmx::power5(2));
    EXPECT_EQ(64,  gmx::power6(2));
    EXPECT_EQ(4096, gmx::power12(2));

    EXPECT_EQ(6.25,               gmx::square(2.5));
    EXPECT_EQ(15.625,             gmx::power3(2.5));
    EXPECT_EQ(39.0625,            gmx::power4(2.5));
    EXPECT_EQ(97.65625,           gmx::power5(2.5));
    EXPECT_EQ(244.140625,         gmx::power6(2.5));
    EXPECT_EQ(59604.644775390625, gmx::power12(2.5));
}

} // namespace
