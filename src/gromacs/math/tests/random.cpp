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
 * \brief Tests for GROMACS random engines and distributions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/random.h"

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace
{

class RandomTest : public ::testing::TestWithParam<std::vector<gmx_uint64_t> >
{
};

TEST_P(RandomTest, ThreeFry2x64)
{
    gmx::test::TestReferenceData       data;
    gmx::test::TestReferenceChecker    checker(data.rootChecker());
    const std::vector<gmx_uint64_t>    input = GetParam();
    std::vector<gmx_uint64_t>          result;

    gmx::ThreeFry2x64<0>               rng(input[2], input[3]);
    rng.restart(input[0], input[1]);

    result.push_back(rng());
    result.push_back(rng());

    checker.checkSequence(result.begin(), result.end(), "ThreeFry2x64");
}

TEST_P(RandomTest, ThreeFry2x64Fast)
{
    gmx::test::TestReferenceData       data;
    gmx::test::TestReferenceChecker    checker(data.rootChecker());
    const std::vector<gmx_uint64_t>    input = GetParam();
    std::vector<gmx_uint64_t>          result;

    gmx::ThreeFry2x64Fast<0>           rng(input[2], input[3]);
    rng.restart(input[0], input[1]);

    result.push_back(rng());
    result.push_back(rng());

    checker.checkSequence(result.begin(), result.end(), "ThreeFry2x64Fast");
}

TEST_P(RandomTest, ThreeFry2x64Using40Rounds)
{
    gmx::test::TestReferenceData       data;
    gmx::test::TestReferenceChecker    checker(data.rootChecker());
    const std::vector<gmx_uint64_t>    input = GetParam();
    std::vector<gmx_uint64_t>          result;

    gmx::ThreeFry2x64General<40, 0>    rng(input[2], input[3]);
    rng.restart(input[0], input[1]);

    result.push_back(rng());
    result.push_back(rng());

    checker.checkSequence(result.begin(), result.end(), "ThreeFry2x64Using40Rounds");
}


/*! \brief Constant array of integers with all bits zeroed.
 *
 *  Reference key and counter input data for known answers test.
 *  The 2x64 flavors of ThreeFry64 will use the first four values, while
 *  the 4x64 version uses all eight.
 */
const std::vector<gmx_uint64_t> bitsZero {{
                                              0, 0, 0, 0
                                          }};


/*! \brief Constant array of integers with all bits set to one.
 *
 *  Reference key and counter input data for known answers test.
 *  The 2x64 flavors of ThreeFry64 will use the first four values, while
 *  the 4x64 version uses all eight.
 */
const std::vector<gmx_uint64_t> bitsOne {{
                                             0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
                                             0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL
                                         }};

/*! \brief Constant array of integers with bitpattern from Pi.
 *
 *  Reference key and counter input data for known answers test.
 *  The 2x64 flavors of ThreeFry64 will use the first four values, while
 *  the 4x64 version uses all eight.
 */
const std::vector<gmx_uint64_t> bitsPi {{
                                            0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL,
                                            0xa4093822299f31d0ULL, 0x082efa98ec4e6c89ULL
                                        }};

// Test the known ansers for the ThreeFry random function when the argument
// is (1) all zero, (2) all ones, (3) the bits of pi, for a bunch of different flavors of ThreeFry.
INSTANTIATE_TEST_CASE_P(KnownAnswersTest, RandomTest,
                            ::testing::Values(bitsZero, bitsOne, bitsPi));


// Test the random device call
TEST_F(RandomTest, makeRandomSeed)
{
    // Unlike Sony, we do not use "4" as a constant random value, so the only
    // thing we can check for the random device is that multiple calls to
    // it produce different results.
    // We choose to ignore the 2^-64 probability this will happen by chance;
    // if you execute the unit tests once per second you might have to run them
    // an extra time rougly once per 300 billion years - apologies in advance!

    gmx_uint64_t i0 = makeRandomSeed();
    gmx_uint64_t i1 = makeRandomSeed();

    EXPECT_NE(i0, i1);
}



// ThreeFry2x64 tests
TEST_F(RandomTest, ThreeFry2x64Logical)
{
    gmx::ThreeFry2x64<10> rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<10> rngB(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<10> rngC(123456, gmx::RandomDomain::Other);

    rngB();              // draw just once first, so block is the same, but index has changed
    EXPECT_NE(rngA, rngB);
    rngC(); rngC();      // two draws: next block, but index is the same
    EXPECT_NE(rngA, rngC);
    rngA();
    EXPECT_EQ(rngA, rngB);
    rngA();
    EXPECT_EQ(rngA, rngC);
}

TEST_F(RandomTest, ThreeFry2x64InternalCounterSequence)
{
    gmx::test::TestReferenceData     data;
    gmx::test::TestReferenceChecker  checker(data.rootChecker());

    // 66 bits of internal counter means the first four increments (giving 2*4=8 results)
    // correspond to incrementing word 0, and then we should carry over to word 1.
    gmx::ThreeFry2x64<66>        rngA(123456, gmx::RandomDomain::Other);
    std::vector<gmx_uint64_t>    result;

    for (int i = 0; i < 16; i++)
    {
        result.push_back(rngA());
    }
    checker.checkSequence(result.begin(), result.end(), "ThreeFry2x64InternalCounterSequence");

    // Make sure nothing goes wrong with the internal counter sequence when we use a full 64-bit word
    gmx::ThreeFry2x64<64>        rngB(123456, gmx::RandomDomain::Other);
    for (int i = 0; i < 16; i++)
    {
        rngB();
    }

    // Use every single bit for the internal counter
    gmx::ThreeFry2x64<128>        rngC(123456, gmx::RandomDomain::Other);
    for (int i = 0; i < 16; i++)
    {
        rngC();
    }
}

TEST_F(RandomTest, ThreeFry2x64Reseed)
{
    gmx::ThreeFry2x64<10> rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<10> rngB;

    EXPECT_NE(rngA, rngB);
    rngB.seed(123456, gmx::RandomDomain::Other);
    EXPECT_EQ(rngA, rngB);
    rngB();                                      // internal counter increments
    rngB.seed(123456, gmx::RandomDomain::Other); // reseeding should reset random stream too
    EXPECT_EQ(rngA, rngB);
}

TEST_F(RandomTest, ThreeFry2x64Discard)
{
    gmx::ThreeFry2x64<10> rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<10> rngB(123456, gmx::RandomDomain::Other);

    for (int i = 0; i < 9; i++)
    {
        rngA();
    }
    rngB.discard(9);
    EXPECT_EQ(rngA, rngB);
}


TEST_F(RandomTest, ThreeFry2x64InvalidCounter)
{
    gmx::ThreeFry2x64<10> rngA(123456, gmx::RandomDomain::Other);

    // Highest 10 bits of counter reserved for the internal counter.
    EXPECT_THROW_GMX(rngA.restart(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF), gmx::InternalError);
}

TEST_F(RandomTest, ThreeFry2x64ExhaustInternalCounter)
{
    gmx::ThreeFry2x64<2> rngA(123456, gmx::RandomDomain::Other);

    // 2 bits for internal counter and 2 64-results per counter means 8 results are fine
    for (int i = 0; i < 8; i++)
    {
        rngA();
    }
    // ... but the 9th time we have exhausted the internal counter space.
    EXPECT_THROW_GMX(rngA(), gmx::InternalError);
}


TEST_F(RandomTest, TabulatedNormalDistribution14)
{
    gmx::test::TestReferenceData         data;
    gmx::test::TestReferenceChecker      checker(data.rootChecker());

    gmx::ThreeFry2x64<2>                 rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<>   dist(2.0, 5.0); // Use default 14-bit resolution
    std::vector<float>                   result;

    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "TabulatedNormalDistribution14");
}

TEST_F(RandomTest, TabulatedNormalDistribution16)
{
    gmx::test::TestReferenceData                  data;
    gmx::test::TestReferenceChecker               checker(data.rootChecker());

    gmx::ThreeFry2x64<2>                          rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<float, 16>   dist(2.0, 5.0); // Use larger 16-bit table
    std::vector<float>                            result;

    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "TabulatedNormalDistribution10");
}

TEST_F(RandomTest, TabulatedNormalDistributionDouble14)
{
    gmx::test::TestReferenceData                  data;
    gmx::test::TestReferenceChecker               checker(data.rootChecker());

    gmx::ThreeFry2x64<2>                          rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<double>      dist(2.0, 5.0);
    std::vector<double>                           result;

    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "TabulatedNormalDistribution10");
}

TEST_F(RandomTest, TabulatedNormalDistributionLogical)
{
    gmx::ThreeFry2x64<2>                 rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<>   distA(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>   distB(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>   distC(3.0, 5.0);
    gmx::TabulatedNormalDistribution<>   distD(2.0, 4.0);

    EXPECT_EQ(distA, distB);
    EXPECT_NE(distA, distC);
    EXPECT_NE(distA, distD);
}


TEST_F(RandomTest, TabulatedNormalDistributionReset)
{
    gmx::ThreeFry2x64<2>                                      rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<>                        distA(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>                        distB(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>::result_type           valA, valB;

    valA = distA(rng);

    distB(rng);
    rng.restart();
    distB.reset();

    valB = distB(rng);

    EXPECT_EQ(valA, valB);
}

TEST_F(RandomTest, TabulatedNormalDistributionAltParam)
{
    gmx::ThreeFry2x64<2>                            rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<2>                            rngB(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<>              distA(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>              distB;
    gmx::TabulatedNormalDistribution<>::param_type  paramA(2.0, 5.0);

    EXPECT_NE(distA(rngA), distB(rngB));
    rngA.restart();
    rngB.restart();
    distA.reset();
    distB.reset();
    EXPECT_EQ(distA(rngA), distB(rngB, paramA));
}



TEST_F(RandomTest, GammaDistribution)
{
    gmx::test::TestReferenceData       data;
    gmx::test::TestReferenceChecker    checker(data.rootChecker());

    gmx::ThreeFry2x64<8>               rng(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real>       dist(2.0, 5.0);
    std::vector<real>                  result;

    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "GammaDistribution");
}

TEST_F(RandomTest, GammaDistributionLogical)
{
    gmx::ThreeFry2x64<8>           rng(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real>   distA(2.0, 5.0);
    gmx::GammaDistribution<real>   distB(2.0, 5.0);
    gmx::GammaDistribution<real>   distC(3.0, 5.0);
    gmx::GammaDistribution<real>   distD(2.0, 4.0);

    EXPECT_EQ(distA, distB);
    EXPECT_NE(distA, distC);
    EXPECT_NE(distA, distD);
}


TEST_F(RandomTest, GammaDistributionReset)
{
    gmx::ThreeFry2x64<8>                                rng(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real>                        distA(2.0, 5.0);
    gmx::GammaDistribution<real>                        distB(2.0, 5.0);
    gmx::GammaDistribution<>::result_type               valA, valB;

    valA = distA(rng);

    distB(rng);
    rng.restart();
    distB.reset();

    valB = distB(rng);

    EXPECT_EQ(valA, valB);
}

TEST_F(RandomTest, GammaDistributionAltParam)
{
    gmx::ThreeFry2x64<8>                      rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<8>                      rngB(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real>              distA(2.0, 5.0);
    gmx::GammaDistribution<real>              distB; // default parameters
    gmx::GammaDistribution<real>::param_type  paramA(2.0, 5.0);

    EXPECT_NE(distA(rngA), distB(rngB));
    rngA.restart();
    rngB.restart();
    distA.reset();
    distB.reset();
    EXPECT_EQ(distA(rngA), distB(rngB, paramA));
}



}      // namespace anonymous

}      // namespace gmx
