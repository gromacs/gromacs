/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Tests for SYCL-kernel functionality
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/sycl_kernel_utils.h"

#include <numeric>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/template_mp.h"

#include "testutils/test_hardware_environment.h"

namespace gmx
{
namespace test
{
namespace
{

template<bool useConditionalStaticLocalStorage>
auto buildKernel(sycl::handler& cgh, const int* gm_input, int* gm_output)
{
    // Organize static local storage when using sycl::local_accessor
    using Buffer            = StaticLocalStorage<int, 1>;
    using ConditionalBuffer = StaticLocalStorage<int, 1, useConditionalStaticLocalStorage>;
    // These declarations must be made in the host code
    auto sm_bufferHostStorage            = Buffer::makeHostStorage(cgh);
    auto sm_conditionalBufferHostStorage = ConditionalBuffer::makeHostStorage(cgh);

    return [=](sycl::nd_item<1> itemNdIdx)
    {
        int itemIdx = itemNdIdx.get_global_linear_id();
        // These declarations work in the device kernel.
        typename Buffer::DeviceStorage            sm_bufferDeviceStorage;
        typename ConditionalBuffer::DeviceStorage sm_conditionalBufferDeviceStorage;
        // Extract the valid pointers to local storage
        sycl::local_ptr<int> sm_buffer = Buffer::get_pointer(sm_bufferHostStorage, sm_bufferDeviceStorage);
        sycl::local_ptr<int> sm_conditionalBuffer = ConditionalBuffer::get_pointer(
                sm_conditionalBufferHostStorage, sm_conditionalBufferDeviceStorage);
        sm_buffer[0] = gm_input[itemIdx];
        if constexpr (useConditionalStaticLocalStorage)
        {
            sm_conditionalBuffer[0] = sm_buffer[0];
            gm_output[itemIdx]      = sm_conditionalBuffer[0];
        }
        else
        {
            gm_output[itemIdx] = sm_buffer[0];
        }
    };
}

//! Test fixture parameterized on whether to test conditional static local storage as well
struct StaticLocalStorageTest : public ::testing::TestWithParam<bool>
{
};

/*! \brief Help GoogleTest name our test cases */
std::string nameOfTest(const testing::TestParamInfo<bool>& info)
{
    std::string testName = info.param ? "use_conditional_storage" : "no_conditional_storage";

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");
    return testName;
}

TEST_P(StaticLocalStorageTest, Works)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        const DeviceContext& deviceContext = testDevice->deviceContext();
        const DeviceStream&  deviceStream  = testDevice->deviceStream();
        deviceContext.activate();

        // Prepare inputs
        const int       numThreads  = 32;
        const int       lowestValue = 3;
        HostVector<int> input(numThreads, { PinningPolicy::PinnedIfSupported });
        // Give each thread a unique value to handle
        std::iota(input.begin(), input.end(), lowestValue);
        DeviceBuffer<int> d_input;
        allocateDeviceBuffer(&d_input, input.size(), deviceContext);
        copyToDeviceBuffer(
                &d_input, input.data(), 0, input.size(), deviceStream, GpuApiCallBehavior::Sync, nullptr);
        const int* gm_input = d_input.get_pointer();

        // Prepare outputs, using initial sentinel value of -1
        HostVector<int>   output(input.size(), -1, { PinningPolicy::PinnedIfSupported });
        DeviceBuffer<int> d_output;
        allocateDeviceBuffer(&d_output, output.size(), deviceContext);
        int* gm_output = d_output.get_pointer();

        // Submit device kernel to copy from input to static local storage to output
        dispatchTemplatedFunction(
                [&](auto useConditionalStaticLocalStorage)
                {
                    syclSubmitWithoutEvent(
                            deviceStream.stream(),
                            [&](sycl::handler& cgh)
                            {
                                auto kernel = buildKernel<useConditionalStaticLocalStorage>(
                                        cgh, gm_input, gm_output);
                                cgh.parallel_for(sycl::nd_range<1>(numThreads, numThreads), kernel);
                            });
                },
                GetParam());
        deviceStream.synchronize();

        // Get result from device
        copyFromDeviceBuffer(
                output.data(), &d_output, 0, output.size(), deviceStream, GpuApiCallBehavior::Sync, nullptr);

        // Check result
        std::vector<int> expectedValues(numThreads);
        std::iota(expectedValues.begin(), expectedValues.end(), int(lowestValue));
        EXPECT_THAT(output, testing::Pointwise(testing::Eq(), expectedValues));

        // Clean up
        freeDeviceBuffer(&d_input);
        freeDeviceBuffer(&d_output);
    }
}

INSTANTIATE_TEST_SUITE_P(WithParameters, StaticLocalStorageTest, ::testing::Values(false, true), nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
