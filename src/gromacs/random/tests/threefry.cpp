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
 * \brief Tests for the ThreeFry random engine
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \ingroup module_random
 */
#include "gmxpre.h"

#include "gromacs/random/threefry.h"

#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/random/seed.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "threefrytestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Add runners to the provided list. GPU runners are added if they are supported.
 *
 * \tparam HostTestRunner The type of Threefry test runner to use for the host tests.
 * \tparam DeviceTestRunner The type of ThreeFry test runner to use for the device tests.
 * \param[in,out] runners The list of runners to populate. Should be a vector of unique pointers.
 * \param[in] key0 The first 64 bits of the 128 bit seed.
 * \param[in] key1 The last 64 bits of the 128 bit seed.
 */
template<typename HostTestRunner, typename DeviceTestRunner, typename T>
void generateListOfRunners(T& runners, uint64_t key0, uint64_t key1)
{
    // Add runners for the CPU version
    runners.emplace_back(std::make_unique<HostTestRunner>(key0, key1));
    // If supported, add runners for the GPU version for each available GPU
    constexpr bool addGpuRunners = GPU_THREEFRY_SUPPORTED;
    if (addGpuRunners)
    {
        for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
        {
            runners.emplace_back(std::make_unique<DeviceTestRunner>(*testDevice, key0, key1));
        }
    }
}

/*! \brief Loop over a list of test runners, e.g., on CPU and GPUs, for the same test.
 *         Compare the results to the common reference data.
 * \param[in] runners The list of test runner to run
 * \param[in] testName The name of the test
 * \param[in] input The test input
 */
template<typename T>
void runAllTestRunners(const T&                                      runners,
                       [[maybe_unused]] const std::string&           testName,
                       [[maybe_unused]] const std::vector<uint64_t>& input)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    for (const auto& runner : runners)
    {
        std::string testDescription = formatString("Testing %s on %s with input: %lu %lu %lu %lu",
                                                   testName.c_str(),
                                                   runner->hardwareDescription().c_str(),
                                                   input[0],
                                                   input[1],
                                                   input[2],
                                                   input[3]);
        SCOPED_TRACE(testDescription);
        std::vector<uint64_t> result;

        runner->restartRng(input[0], input[1]);

        result.push_back(runner->nextRng());
        result.push_back(runner->nextRng());

        checker.checkSequence(result.begin(), result.end(), testName.c_str());
    }
}

class ThreeFry2x64Test : public ::testing::TestWithParam<std::vector<uint64_t>>
{
};

TEST_P(ThreeFry2x64Test, Default)
{
    // Construct the list of runners
    std::vector<std::unique_ptr<IThreeFryGeneralTestRunner<20, 0>>> runners;

    const std::vector<uint64_t> input = GetParam();
    // Add CPU and GPU runners to the list
    generateListOfRunners<ThreeFryHostTestRunner<0>, ThreeFryDeviceTestRunner<0>>(
            runners, input[2], input[3]);

    runAllTestRunners(runners, "ThreeFry2x64", input);
}

TEST_P(ThreeFry2x64Test, Fast)
{
    // Construct the list of runners
    std::vector<std::unique_ptr<IThreeFryGeneralTestRunner<13, 0>>> runners;

    const std::vector<uint64_t> input = GetParam();
    // Add CPU and GPU runners to the list
    generateListOfRunners<ThreeFryFastHostTestRunner<0>, ThreeFryFastDeviceTestRunner<0>>(
            runners, input[2], input[3]);

    runAllTestRunners(runners, "ThreeFry2x64Fast", input);
}

TEST_P(ThreeFry2x64Test, Using40Rounds)
{
    // Construct the list of runners
    std::vector<std::unique_ptr<IThreeFryGeneralTestRunner<40, 0>>> runners;

    const std::vector<uint64_t> input = GetParam();
    // Add CPU and GPU runners to the list
    generateListOfRunners<ThreeFryGeneralHostTestRunner<40, 0>, ThreeFryGeneralDeviceTestRunner<40, 0>>(
            runners, input[2], input[3]);

    runAllTestRunners(runners, "ThreeFry2x64Using40Rounds", input);
}


/*! \brief Constant array of integers with all bits zeroed.
 *
 *  Reference key and counter input data for known answers test.
 *  The 2x64 flavors of ThreeFry64 will use the first four values, while
 *  the 4x64 version uses all eight.
 */
const std::vector<uint64_t> bitsZero{ { 0, 0, 0, 0 } };


/*! \brief Constant array of integers with all bits set to one.
 *
 *  Reference key and counter input data for known answers test.
 *  The 2x64 flavors of ThreeFry64 will use the first four values, while
 *  the 4x64 version uses all eight.
 */
const std::vector<uint64_t> bitsOne{
    { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL }
};

/*! \brief Constant array of integers with bitpattern from Pi.
 *
 *  Reference key and counter input data for known answers test.
 *  The 2x64 flavors of ThreeFry64 will use the first four values, while
 *  the 4x64 version uses all eight.
 */
const std::vector<uint64_t> bitsPi{
    { 0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL, 0xa4093822299f31d0ULL, 0x082efa98ec4e6c89ULL }
};

// Test the known answers for the ThreeFry random function when the argument
// is (1) all zero, (2) all ones, (3) the bits of pi, for a bunch of different flavors of ThreeFry.
// This suite of tests are run on both CPU and GPU (if available).
INSTANTIATE_TEST_SUITE_P(KnownAnswersTest, ThreeFry2x64Test, ::testing::Values(bitsZero, bitsOne, bitsPi));


// ThreeFry2x64 tests
TEST_F(ThreeFry2x64Test, Logical)
{
    gmx::ThreeFry2x64<10> rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<10> rngB(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<10> rngC(123456, gmx::RandomDomain::Other);

    rngB(); // draw just once first, so block is the same, but index has changed
    EXPECT_NE(rngA, rngB);
    rngC();
    rngC(); // two draws: next block, but index is the same
    EXPECT_NE(rngA, rngC);
    rngA();
    EXPECT_EQ(rngA, rngB);
    rngA();
    EXPECT_EQ(rngA, rngC);
}

TEST_F(ThreeFry2x64Test, InternalCounterSequence)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // 66 bits of internal counter means the first four increments (giving 2*4=8 results)
    // correspond to incrementing word 0, and then we should carry over to word 1.
    gmx::ThreeFry2x64<66> rngA(123456, gmx::RandomDomain::Other);
    std::vector<uint64_t> result;

    result.reserve(16);
    for (int i = 0; i < 16; i++)
    {
        result.push_back(rngA());
    }
    checker.checkSequence(result.begin(), result.end(), "ThreeFry2x64InternalCounterSequence");

    // Make sure nothing goes wrong with the internal counter sequence when we use a full 64-bit word
    gmx::ThreeFry2x64<64> rngB(123456, gmx::RandomDomain::Other);
    for (int i = 0; i < 16; i++)
    {
        rngB();
    }

    // Use every single bit for the internal counter
    gmx::ThreeFry2x64<128> rngC(123456, gmx::RandomDomain::Other);
    for (int i = 0; i < 16; i++)
    {
        rngC();
    }
}

TEST_F(ThreeFry2x64Test, Reseed)
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

TEST_F(ThreeFry2x64Test, Discard)
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


TEST_F(ThreeFry2x64Test, InvalidCounter)
{
    gmx::ThreeFry2x64<10> rngA(123456, gmx::RandomDomain::Other);

    // Highest 10 bits of counter reserved for the internal counter.
    EXPECT_THROW_GMX(rngA.restart(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF), gmx::InternalError);
}

TEST_F(ThreeFry2x64Test, ExhaustInternalCounter)
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

} // namespace

} // namespace test
} // namespace gmx
