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
 * \brief Basic tests for working GPU-aware MPI
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "config.h"

#include <numeric>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#if GMX_GPU
#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/hostallocator.h"
#endif

#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/test_hardware_environment.h"

namespace gmx
{

namespace test
{

namespace
{

MessageStringCollector getSkipMessagesIfNecessary()
{
    MessageStringCollector errorReasons;

    errorReasons.appendIf(!GMX_LIB_MPI,
                          "Can only test GPU-aware MPI on a build with an MPI library");
    errorReasons.appendIf(!GMX_GPU, "Can only test GPU-aware MPI in a build with GPU support");
    errorReasons.appendIf(getTestHardwareEnvironment()->getTestDeviceList().empty(),
                          "Can only test GPU-aware MPI on ranks with available GPUs");
    errorReasons.appendIf(getTestHardwareEnvironment()->hwinfo()->minGpuAwareMpiStatus
                                  == GpuAwareMpiStatus::NotSupported,
                          "GPU-aware MPI is not supported on at least one device and/or rank");

    return errorReasons;
}

struct GpuAwareMpiTestParams
{
    //! The test-parameter tuple, must match the layout of the rest of the class
    using TupleT = std::tuple<int>;
    //! Number of values present in the MPI message
    int numValuesInMessage;
    GpuAwareMpiTestParams(TupleT t) : numValuesInMessage(std::get<0>(t)) {}
};

/*! \brief Help GoogleTest name our test cases
 *
 * If changes are needed here, consider making matching changes in
 * makeRefDataFileName(). */
std::string nameOfTest(const testing::TestParamInfo<GpuAwareMpiTestParams>& info)
{
    std::string testName = formatString("%d_values_in_message", info.param.numValuesInMessage);

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

//! Test fixture for GPU-aware MPI
class GpuAwareMpiTest : public ::testing::TestWithParam<GpuAwareMpiTestParams>
{
public:
    MpiComm   mpiComm_{ MPI_COMM_WORLD };
    const int unpinnedHostToDeviceTag_ = 0;
    const int pinnedHostToDeviceTag_   = 1;
    const int deviceToUnpinnedHostTag_ = 2;
    const int deviceToPinnedHostTag_   = 3;
    const int deviceToDeviceTag_       = 4;
    const int offset_                  = 0;
    const int specialValue_            = 42;

#if GMX_GPU
    void sendFromDeviceBuffer(const int partnerRank, const int mpiTag)
    {
        const TestDevice* testDevice = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
        const DeviceContext& deviceContext = testDevice->deviceContext();
        const DeviceStream&  deviceStream  = testDevice->deviceStream();
        deviceContext.activate();

        const auto [numValuesInMessage] = GetParam();

        std::vector<int> valuesToSend(numValuesInMessage);
        std::iota(valuesToSend.begin(), valuesToSend.end(), specialValue_ + mpiComm_.rank());

        DeviceBuffer<int> sendBuffer;
        allocateDeviceBuffer(&sendBuffer, numValuesInMessage, deviceContext);
        // Copy H2D
        copyToDeviceBuffer(&sendBuffer,
                           valuesToSend.data(),
                           offset_,
                           numValuesInMessage,
                           deviceStream,
                           GpuApiCallBehavior::Sync,
                           nullptr);

        // Send to partner rank
        MPI_Send(asMpiPointer(sendBuffer), numValuesInMessage, MPI_INT, partnerRank, mpiTag, mpiComm_.comm());

        freeDeviceBuffer(&sendBuffer);
    }

    void receiveToDeviceBuffer(const int partnerRank, const int mpiTag)
    {
        const TestDevice* testDevice = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
        const DeviceContext& deviceContext = testDevice->deviceContext();
        const DeviceStream&  deviceStream  = testDevice->deviceStream();
        deviceContext.activate();

        const auto [numValuesInMessage] = GetParam();

        DeviceBuffer<int> receiveBuffer;
        allocateDeviceBuffer(&receiveBuffer, numValuesInMessage, deviceContext);
        HostVector<int> valuesReceived(numValuesInMessage, -1, { PinningPolicy::PinnedIfSupported });

        // Post receive from partner rank
        MPI_Recv(asMpiPointer(receiveBuffer),
                 numValuesInMessage,
                 MPI_INT,
                 partnerRank,
                 mpiTag,
                 mpiComm_.comm(),
                 MPI_STATUS_IGNORE);

        // Copy D2H
        copyFromDeviceBuffer(valuesReceived.data(),
                             &receiveBuffer,
                             offset_,
                             numValuesInMessage,
                             deviceStream,
                             GpuApiCallBehavior::Sync,
                             nullptr);
        // Check result
        std::vector<int> expectedValues(numValuesInMessage);
        std::iota(expectedValues.begin(), expectedValues.end(), specialValue_ + partnerRank);
        EXPECT_THAT(valuesReceived, testing::Pointwise(testing::Eq(), expectedValues));

        freeDeviceBuffer(&receiveBuffer);
    }
#endif
};

TEST_P(GpuAwareMpiTest, SendFromUnpinnedHostBufferToDeviceBuffer)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    MessageStringCollector skipReasons = getSkipMessagesIfNecessary();
    if (!skipReasons.isEmpty())
    {
        GTEST_SKIP() << skipReasons.toString();
    }

#if GMX_GPU
    const TestDevice*    testDevice    = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
    const DeviceContext& deviceContext = testDevice->deviceContext();
    deviceContext.activate();

    const int partnerRank           = mpiComm_.size() - 1 - mpiComm_.rank();
    const auto [numValuesInMessage] = GetParam();

    const int mpiTag = unpinnedHostToDeviceTag_;
    if (mpiComm_.rank() == 0)
    {
        // This is used to send x from PP ranks to PME ranks after pair search
        std::vector<int> sendBuffer(numValuesInMessage);
        std::iota(sendBuffer.begin(), sendBuffer.end(), specialValue_ + mpiComm_.rank());

        // Send to partner rank
        MPI_Send(sendBuffer.data(), numValuesInMessage, MPI_INT, partnerRank, mpiTag, mpiComm_.comm());
    }
    else
    {
        receiveToDeviceBuffer(partnerRank, mpiTag);
    }
#endif
}

TEST_P(GpuAwareMpiTest, SendFromPinnedHostBufferToDeviceBuffer)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    MessageStringCollector skipReasons = getSkipMessagesIfNecessary();
    if (!skipReasons.isEmpty())
    {
        GTEST_SKIP() << skipReasons.toString();
    }

#if GMX_GPU
    const TestDevice*    testDevice    = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
    const DeviceContext& deviceContext = testDevice->deviceContext();
    deviceContext.activate();

    const int partnerRank           = mpiComm_.size() - 1 - mpiComm_.rank();
    const auto [numValuesInMessage] = GetParam();

    const int mpiTag = pinnedHostToDeviceTag_;
    if (mpiComm_.rank() == 0)
    {
        // This is used to send x from PP ranks to PME ranks after pair search
        HostVector<int> sendBuffer(numValuesInMessage, { PinningPolicy::PinnedIfSupported });
        std::iota(sendBuffer.begin(), sendBuffer.end(), specialValue_ + mpiComm_.rank());

        // Send to partner rank
        MPI_Send(sendBuffer.data(), numValuesInMessage, MPI_INT, partnerRank, mpiTag, mpiComm_.comm());
    }
    else
    {
        receiveToDeviceBuffer(partnerRank, mpiTag);
    }
#endif
}

TEST_P(GpuAwareMpiTest, SendFromDeviceBufferToUnpinnedHostBuffer)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    MessageStringCollector skipReasons = getSkipMessagesIfNecessary();
    if (!skipReasons.isEmpty())
    {
        GTEST_SKIP() << skipReasons.toString();
    }

#if GMX_GPU
    const TestDevice*    testDevice    = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
    const DeviceContext& deviceContext = testDevice->deviceContext();
    deviceContext.activate();

    const int partnerRank           = mpiComm_.size() - 1 - mpiComm_.rank();
    const auto [numValuesInMessage] = GetParam();

    const int mpiTag = deviceToUnpinnedHostTag_;
    if (mpiComm_.rank() == 0)
    {
        // This is used by settle to return the virial
        sendFromDeviceBuffer(partnerRank, mpiTag);
    }
    else
    {
        std::vector<int> receiveBuffer(numValuesInMessage, -1);

        // Post receive from partner rank
        MPI_Recv(receiveBuffer.data(), numValuesInMessage, MPI_INT, partnerRank, mpiTag, mpiComm_.comm(), MPI_STATUS_IGNORE);

        // Check result
        std::vector<int> expectedValues(numValuesInMessage);
        std::iota(expectedValues.begin(), expectedValues.end(), specialValue_ + partnerRank);
        EXPECT_THAT(receiveBuffer, testing::Pointwise(testing::Eq(), expectedValues));
    }
#endif
}

TEST_P(GpuAwareMpiTest, DISABLED_SendFromDeviceBufferToPinnedHostBuffer)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    MessageStringCollector skipReasons = getSkipMessagesIfNecessary();
    if (!skipReasons.isEmpty())
    {
        GTEST_SKIP() << skipReasons.toString();
    }

#if GMX_GPU
    const TestDevice*    testDevice    = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
    const DeviceContext& deviceContext = testDevice->deviceContext();
    deviceContext.activate();

    const int partnerRank           = mpiComm_.size() - 1 - mpiComm_.rank();
    const auto [numValuesInMessage] = GetParam();

    const int mpiTag = deviceToPinnedHostTag_;
    if (mpiComm_.rank() == 0)
    {
        // This is used by PME ranks to return the force on
        // energy+virial steps when not using staged
        // communication.
        sendFromDeviceBuffer(partnerRank, mpiTag);
    }
    else
    {
        HostVector<int> receiveBuffer(numValuesInMessage, -1, { PinningPolicy::PinnedIfSupported });

        // Post receive from partner rank
        MPI_Recv(receiveBuffer.data(), numValuesInMessage, MPI_INT, partnerRank, mpiTag, mpiComm_.comm(), MPI_STATUS_IGNORE);

        // Check result
        std::vector<int> expectedValues(numValuesInMessage);
        std::iota(expectedValues.begin(), expectedValues.end(), specialValue_ + partnerRank);
        EXPECT_THAT(receiveBuffer, testing::Pointwise(testing::Eq(), expectedValues));
    }
#endif
}

TEST_P(GpuAwareMpiTest, SendFromDeviceBufferToDeviceBuffer)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    MessageStringCollector skipReasons = getSkipMessagesIfNecessary();
    if (!skipReasons.isEmpty())
    {
        GTEST_SKIP() << skipReasons.toString();
    }

#if GMX_GPU
    const TestDevice*    testDevice    = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
    const DeviceContext& deviceContext = testDevice->deviceContext();
    deviceContext.activate();

    const int partnerRank = mpiComm_.size() - 1 - mpiComm_.rank();

    const int mpiTag = deviceToDeviceTag_;
    if (mpiComm_.rank() == 0)
    {
        // This is used in PP-PME communication on most steps
        sendFromDeviceBuffer(partnerRank, mpiTag);
    }
    else
    {
        receiveToDeviceBuffer(partnerRank, mpiTag);
    }
#endif
}

// Test the pattern used in GPU-aware halo exchange
TEST_P(GpuAwareMpiTest, IrecvSendPair)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    MessageStringCollector skipReasons = getSkipMessagesIfNecessary();
    if (!skipReasons.isEmpty())
    {
        GTEST_SKIP() << skipReasons.toString();
    }

#if GMX_GPU
    const TestDevice*    testDevice    = getTestHardwareEnvironment()->getTestDeviceList()[0].get();
    const DeviceContext& deviceContext = testDevice->deviceContext();
    const DeviceStream&  deviceStream  = testDevice->deviceStream();
    deviceContext.activate();

    const int partnerRank           = mpiComm_.size() - 1 - mpiComm_.rank();
    const auto [numValuesInMessage] = GetParam();
    HostVector<int> valuesToSend(numValuesInMessage, { PinningPolicy::PinnedIfSupported });
    HostVector<int> valuesReceived(numValuesInMessage, -1, { PinningPolicy::PinnedIfSupported });

    std::iota(valuesToSend.begin(), valuesToSend.end(), specialValue_ + mpiComm_.rank());
    std::vector<int> expectedValues(numValuesInMessage);
    std::iota(expectedValues.begin(), expectedValues.end(), specialValue_ + partnerRank);

    DeviceBuffer<int> sendBuffer, receiveBuffer;
    allocateDeviceBuffer(&receiveBuffer, numValuesInMessage, deviceContext);

    // Post initial non-blocking receive from partner rank
    MPI_Request request;
    MPI_Irecv(asMpiPointer(receiveBuffer),
              numValuesInMessage,
              MPI_INT,
              partnerRank,
              deviceToDeviceTag_,
              mpiComm_.comm(),
              &request);

    // Copy H2D
    allocateDeviceBuffer(&sendBuffer, numValuesInMessage, deviceContext);
    copyToDeviceBuffer(&sendBuffer,
                       valuesToSend.data(),
                       offset_,
                       numValuesInMessage,
                       deviceStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);

    // Send to partner rank
    MPI_Send(asMpiPointer(sendBuffer), numValuesInMessage, MPI_INT, partnerRank, deviceToDeviceTag_, mpiComm_.comm());
    // Wait for non-blocking receive
    MPI_Wait(&request, MPI_STATUS_IGNORE);

    // Copy D2H
    copyFromDeviceBuffer(valuesReceived.data(),
                         &receiveBuffer,
                         offset_,
                         numValuesInMessage,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);
    // Check result
    EXPECT_THAT(valuesReceived, testing::Pointwise(testing::Eq(), expectedValues));

    freeDeviceBuffer(&sendBuffer);
    freeDeviceBuffer(&receiveBuffer);
#endif
}

INSTANTIATE_TEST_SUITE_P(WorksOnParameters,
                         GpuAwareMpiTest,
                         ::testing::ConvertGenerator<GpuAwareMpiTestParams::TupleT>(::testing::Combine(
                                 // Ensure that zero-, small-, and medium-sized messages work
                                 ::testing::Values(0, 1, 12000))),
                         nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
