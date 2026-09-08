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

#include "config.h"

#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/random/seed.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "threefrytestinputs.h"

#if GMX_GPU && !GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    if GMX_GPU_SYCL
#        include "gromacs/gpu_utils/gmxsycl.h"
#    endif
#endif

namespace gmx
{
namespace test
{
namespace
{

#if GMX_GPU && !GMX_GPU_OPENCL

#    if GMX_GPU_CUDA || GMX_GPU_HIP

template<unsigned int rounds>
GMX_KERNEL_ATTRIBUTE void setupDeviceRngKeys(gmx::ThreeFry2x64General<rounds, 0>* d_rng,
                                             uint64_t                             key0,
                                             uint64_t                             key1)
{
    new (d_rng) gmx::ThreeFry2x64General<rounds, 0>(key0, key1);
}

template<unsigned int rounds>
GMX_KERNEL_ATTRIBUTE void restartDeviceRng(gmx::ThreeFry2x64General<rounds, 0>* d_rng, uint64_t ctr0, uint64_t ctr1)
{
    d_rng->restart(ctr0, ctr1);
}

template<unsigned int rounds>
GMX_KERNEL_ATTRIBUTE void nextDeviceRng(gmx::ThreeFry2x64General<rounds, 0>* d_rng,
                                        uint64_t* __restrict__ gm_result)
{
    *gm_result = (*d_rng)();
}

#    elif GMX_GPU_SYCL

template<unsigned int rounds>
class ThreeFrySetupKernel;
template<unsigned int rounds>
class ThreeFryRestartKernel;
template<unsigned int rounds>
class ThreeFryNextKernel;

#    endif

/*! \brief Run a known-answer ThreeFry test on GPU. */
template<unsigned int rounds>
void runKnownAnswerGpu(const DeviceContext&   deviceContext,
                       const DeviceStream&    deviceStream,
                       uint64_t               ctr0,
                       uint64_t               ctr1,
                       uint64_t               key0,
                       uint64_t               key1,
                       std::vector<uint64_t>& result)
{
    DeviceBuffer<gmx::ThreeFry2x64General<rounds, 0>> rng;
    DeviceBuffer<uint64_t>                            d_result;
    allocateDeviceBuffer(&rng, 1, deviceContext);
    allocateDeviceBuffer(&d_result, 1, deviceContext);

#    if GMX_GPU_CUDA || GMX_GPU_HIP
    {
        KernelLaunchConfig kernelLaunchConfig;

        auto       kernelPtr = setupDeviceRngKeys<rounds>;
        const auto kernelArgs =
                prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &rng, &key0, &key1);
        launchGpuKernel(kernelPtr,
                        kernelLaunchConfig,
                        deviceStream,
                        nullptr,
                        "test_threefry_setup_device_rng_keys_kernel",
                        kernelArgs);
    }
    {
        KernelLaunchConfig kernelLaunchConfig;

        auto       kernelPtr = restartDeviceRng<rounds>;
        const auto kernelArgs =
                prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &rng, &ctr0, &ctr1);
        launchGpuKernel(
                kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "test_threefry_restart_kernel", kernelArgs);
    }
    for (int i = 0; i < 2; ++i)
    {
        KernelLaunchConfig kernelLaunchConfig;

        auto kernelPtr = nextDeviceRng<rounds>;
        const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &rng, &d_result);
        launchGpuKernel(
                kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "test_threefry_next_kernel", kernelArgs);

        uint64_t h_result;
        copyFromDeviceBuffer(&h_result, &d_result, 0, 1, deviceStream, GpuApiCallBehavior::Sync, nullptr);
        result.push_back(h_result);
    }
#    elif GMX_GPU_SYCL
    using Rng   = gmx::ThreeFry2x64General<rounds, 0>;
    Rng* gm_rng = rng.get_pointer();
    deviceStream.stream().submit(
            [&](sycl::handler& cgh) {
                cgh.single_task<ThreeFrySetupKernel<rounds>>([=]() { new (gm_rng) Rng(key0, key1); });
            });
    deviceStream.stream().submit(
            [&](sycl::handler& cgh) {
                cgh.single_task<ThreeFryRestartKernel<rounds>>([=]() { gm_rng->restart(ctr0, ctr1); });
            });
    for (int i = 0; i < 2; ++i)
    {
        uint64_t* gm_dest = d_result.get_pointer();
        deviceStream.stream().submit(
                [&](sycl::handler& cgh)
                { cgh.single_task<ThreeFryNextKernel<rounds>>([=]() { *gm_dest = (*gm_rng)(); }); });

        uint64_t h_result;
        copyFromDeviceBuffer(&h_result, &d_result, 0, 1, deviceStream, GpuApiCallBehavior::Sync, nullptr);
        result.push_back(h_result);
    }
#    endif

    freeDeviceBuffer(&rng);
    freeDeviceBuffer(&d_result);
}

#endif // GMX_GPU && !GMX_GPU_OPENCL

using ThreeFryKnownAnswersTestHelper =
        HardwareAndExecutionTestHelper<ThreeFryKnownAnswerInput, std::tuple<>>;

static const auto sc_knownAnswerInputFormatters = std::make_tuple(
        useString,
        [](uint64_t ctr0) { return formatString("%zu", ctr0).substr(0, 3); },
        [](uint64_t ctr1) { return formatString("%zu", ctr1).substr(0, 3); },
        [](uint64_t key0) { return formatString("%zu", key0).substr(0, 3); },
        [](uint64_t key1) { return formatString("%zu", key1).substr(0, 3); });

static const NameOfTestFromTuple<ThreeFryKnownAnswersTestHelper::DynamicParameters> sc_knownAnswersTestNamer =
        ThreeFryKnownAnswersTestHelper::testNamer(sc_knownAnswerInputFormatters, std::make_tuple());

/*! \brief Test fixture for ThreeFry known-answer tests on CPU and GPU. */
class ThreeFryKnownAnswersTest : public HardwareTestFixture<ThreeFryKnownAnswersTestHelper>
{
protected:
    ThreeFryKnownAnswersTest() : HardwareTestFixture(sc_knownAnswerInputFormatters) {}

    /*! \brief Run a known-answer test on the selected hardware. */
    template<unsigned int rounds>
    void runKnownAnswerTest(const char* sequenceName);
};

template<unsigned int rounds>
void ThreeFryKnownAnswersTest::runKnownAnswerTest(const char* sequenceName)
{
    TestReferenceChecker&      testChecker         = checker();
    const TestHardwareContext* testHardwareContext = hardwareContext();
    auto [inputName, ctr0, ctr1, key0, key1, _]    = GetParam();

    SCOPED_TRACE(formatString("Testing %sUsing%uRounds on %s with input %s",
                              sequenceName,
                              rounds,
                              hardwareContext()->description().c_str(),
                              inputName.c_str()));

    std::vector<uint64_t> result;

    if (testHardwareContext->isGpuTest())
    {
#if GMX_GPU && !GMX_GPU_OPENCL
        testHardwareContext->activate();
        runKnownAnswerGpu<rounds>(*testHardwareContext->deviceContext(),
                                  *testHardwareContext->deviceStream(),
                                  ctr0,
                                  ctr1,
                                  key0,
                                  key1,
                                  result);
#else
        GMX_THROW(gmx::InternalError("GPU hardware context with no test code"));
#endif
    }
    else
    {
        gmx::ThreeFry2x64General<rounds, 0> rng(key0, key1);
        rng.restart(ctr0, ctr1);
        result.push_back(rng());
        result.push_back(rng());
    }

    testChecker.checkSequence(result.begin(), result.end(), sequenceName);
}

TEST_P(ThreeFryKnownAnswersTest, Default)
{
    runKnownAnswerTest<20>("ThreeFry2x64");
}

TEST_P(ThreeFryKnownAnswersTest, Fast)
{
    runKnownAnswerTest<13>("ThreeFry2x64Fast");
}

TEST_P(ThreeFryKnownAnswersTest, Using40Rounds)
{
    runKnownAnswerTest<40>("ThreeFry2x64Using40Rounds");
}

INSTANTIATE_TEST_SUITE_P(AllHardware,
                         ThreeFryKnownAnswersTest,
                         ::testing::ConvertGenerator(
                                 ::testing::Combine(::testing::ValuesIn(sc_threefryKnownAnswerInputs),
                                                    ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                            GpuConfigurationCapabilities::Threefry))),
                                 flattenTupleWithHardwareContext<ThreeFryKnownAnswerInput>()),
                         sc_knownAnswersTestNamer);


class ThreeFry2x64Test : public ::testing::Test
{
};

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
