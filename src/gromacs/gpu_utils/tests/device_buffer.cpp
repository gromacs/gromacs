/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * \brief Tests for device buffer
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_GPU
#    include <numeric>

#    include <gmock/gmock.h>
#    include <gtest/gtest.h>

#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/hostallocator.h"

#    include "testutils/test_hardware_environment.h"
#    include "testutils/testasserts.h"

namespace gmx
{

template<typename ValueType>
BasicVector<ValueType>& operator++(BasicVector<ValueType>& in)
{
    in[XX]++;
    in[YY]++;
    in[ZZ]++;
    return in;
}

template<typename ValueType>
BasicVector<ValueType>& operator++(BasicVector<ValueType>& in, int /* n */)
{
    BasicVector<ValueType> temp = *in;
    ++*in;
    return temp;
}

template<typename ValueType>
inline bool operator==(const BasicVector<ValueType>& lhs, const BasicVector<ValueType>& rhs)
{
    return lhs[XX] == rhs[XX] && lhs[YY] == rhs[YY] && lhs[ZZ] == rhs[ZZ];
}

namespace test
{

namespace
{

using testing::Eq;
using testing::Pointwise;

//! Test fixture (needed for typed tests)
template<typename T>
class DeviceBufferTest : public ::testing::Test
{
};

using TypeParamList = testing::Types<short, int, float, double, gmx::RVec>;
TYPED_TEST_SUITE(DeviceBufferTest, TypeParamList);

TYPED_TEST(DeviceBufferTest, CanAllocateAndFreeDeviceBuffer)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        const DeviceContext& deviceContext = testDevice->deviceContext();
        setActiveDevice(testDevice->deviceInfo());

        DeviceBuffer<TypeParam> buffer;
        int                     numValues = 123;
        allocateDeviceBuffer(&buffer, numValues, deviceContext);
        freeDeviceBuffer(&buffer);
    }
}

TYPED_TEST(DeviceBufferTest, CanReallocateAndFreeDeviceBuffer)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        const DeviceContext& deviceContext = testDevice->deviceContext();
        setActiveDevice(testDevice->deviceInfo());

        DeviceBuffer<TypeParam> buffer;
        int                     currentNumValues    = 456;
        int                     newNumValues        = 789;
        int                     currentMaxNumValues = 0;
        allocateDeviceBuffer(&buffer, currentNumValues, deviceContext);
        reallocateDeviceBuffer(&buffer, newNumValues, &currentNumValues, &currentMaxNumValues, deviceContext);
        freeDeviceBuffer(&buffer);
    }
}

//! Initial value to fill the buffer of the scalar type
template<typename T>
const T c_initialValue = static_cast<T>(1);

//! Initial value to fill the buffer of the vector type
template<>
const gmx::RVec c_initialValue<gmx::RVec> = { 1, -2, 3 };


TYPED_TEST(DeviceBufferTest, CanCopyToAndFromDevice)
{
    for (auto transferKind : { GpuApiCallBehavior::Sync, GpuApiCallBehavior::Async })
    {
        PinningPolicy pinningPolicy = (transferKind == GpuApiCallBehavior::Async)
                                              ? PinningPolicy::PinnedIfSupported
                                              : PinningPolicy::CannotBePinned;
        for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
        {
            const DeviceContext& deviceContext = testDevice->deviceContext();
            const DeviceStream&  deviceStream  = testDevice->deviceStream();
            setActiveDevice(testDevice->deviceInfo());

            DeviceBuffer<TypeParam> buffer;
            int                     numValues = 123;
            allocateDeviceBuffer(&buffer, numValues, deviceContext);
            HostVector<TypeParam> valuesIn(numValues, { pinningPolicy });
            HostVector<TypeParam> valuesOut(numValues, { pinningPolicy });

            std::iota(valuesIn.begin(), valuesIn.end(), c_initialValue<TypeParam>);

            copyToDeviceBuffer(&buffer, valuesIn.data(), 0, numValues, deviceStream, transferKind, nullptr);
            copyFromDeviceBuffer(
                    valuesOut.data(), &buffer, 0, numValues, deviceStream, transferKind, nullptr);
            if (transferKind == GpuApiCallBehavior::Async)
            {
                deviceStream.synchronize();
            }
            EXPECT_THAT(valuesOut, Pointwise(Eq(), valuesIn))
                    << "Changed after H2D and D2H " << enumValueToString(transferKind) << " copy.";
            freeDeviceBuffer(&buffer);
        }
    }
}

TYPED_TEST(DeviceBufferTest, CanCopyToAndFromDeviceWithOffset)
{
    for (auto transferKind : { GpuApiCallBehavior::Sync, GpuApiCallBehavior::Async })
    {
        PinningPolicy pinningPolicy = (transferKind == GpuApiCallBehavior::Async)
                                              ? PinningPolicy::PinnedIfSupported
                                              : PinningPolicy::CannotBePinned;
        for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
        {
            const DeviceContext& deviceContext = testDevice->deviceContext();
            const DeviceStream&  deviceStream  = testDevice->deviceStream();
            setActiveDevice(testDevice->deviceInfo());

            DeviceBuffer<TypeParam> buffer;
            int                     numValues = 123;
            allocateDeviceBuffer(&buffer, 2 * numValues, deviceContext);
            HostVector<TypeParam> valuesIn(numValues, { pinningPolicy });
            HostVector<TypeParam> valuesOut(2 * numValues, { pinningPolicy });

            std::iota(valuesIn.begin(), valuesIn.end(), c_initialValue<TypeParam>);

            // Fill the buffer with two copies of valuesIn, one after the other.
            copyToDeviceBuffer(&buffer, valuesIn.data(), 0, numValues, deviceStream, transferKind, nullptr);
            copyToDeviceBuffer(
                    &buffer, valuesIn.data(), numValues, numValues, deviceStream, transferKind, nullptr);
            // Wait until GPU is done andd o the same copying on the CPU, so we can test it works correctly.
            if (transferKind == GpuApiCallBehavior::Async)
            {
                deviceStream.synchronize();
            }
            valuesIn.insert(valuesIn.end(), valuesIn.begin(), valuesIn.end());

            copyFromDeviceBuffer(
                    valuesOut.data(), &buffer, 0, 2 * numValues, deviceStream, transferKind, nullptr);
            if (transferKind == GpuApiCallBehavior::Async)
            {
                deviceStream.synchronize();
            }
            EXPECT_THAT(valuesOut, Pointwise(Eq(), valuesIn))
                    << "Changed after H2D and D2H " << enumValueToString(transferKind) << " copy.";

            SCOPED_TRACE("Checking the copy respects the output range");

            // Remove the first element, and push another copy of the last
            // element, so we can check that a copy of all of the data
            // skipping the first element correctly over-writes exactly
            // all but one of the old values.
            valuesIn.erase(valuesIn.begin());
            valuesIn.push_back(valuesIn.back());
            copyFromDeviceBuffer(
                    valuesOut.data(), &buffer, 1, 2 * numValues - 1, deviceStream, transferKind, nullptr);
            if (transferKind == GpuApiCallBehavior::Async)
            {
                deviceStream.synchronize();
            }
            EXPECT_THAT(valuesOut, Pointwise(Eq(), valuesIn))
                    << "Changed after H2D and D2H " << enumValueToString(transferKind) << " copy.";
        }
    }
}

#    if GMX_GPU_CUDA

TYPED_TEST(DeviceBufferTest, CanCopyBetweenDeviceBuffers)
{
    for (auto transferKind : { GpuApiCallBehavior::Sync, GpuApiCallBehavior::Async })
    {
        PinningPolicy pinningPolicy = (transferKind == GpuApiCallBehavior::Async)
                                              ? PinningPolicy::PinnedIfSupported
                                              : PinningPolicy::CannotBePinned;
        for (const auto& testDeviceIn : getTestHardwareEnvironment()->getTestDeviceList())
        {
            for (const auto& testDeviceOut : getTestHardwareEnvironment()->getTestDeviceList())
            {
                int                   numValues = 321;
                HostVector<TypeParam> valuesIn(numValues, { pinningPolicy });
                HostVector<TypeParam> valuesOut(numValues, { pinningPolicy });

                std::iota(valuesIn.begin(), valuesIn.end(), c_initialValue<TypeParam>);

                const DeviceContext& deviceContextIn = testDeviceIn->deviceContext();
                const DeviceStream&  deviceStreamIn  = testDeviceIn->deviceStream();
                setActiveDevice(testDeviceIn->deviceInfo());
                DeviceBuffer<TypeParam> bufferIn;
                allocateDeviceBuffer(&bufferIn, numValues, deviceContextIn);

                const DeviceContext& deviceContextOut = testDeviceOut->deviceContext();
                const DeviceStream&  deviceStreamOut  = testDeviceOut->deviceStream();
                setActiveDevice(testDeviceOut->deviceInfo());
                DeviceBuffer<TypeParam> bufferOut;
                allocateDeviceBuffer(&bufferOut, numValues, deviceContextOut);

                copyToDeviceBuffer(
                        &bufferIn, valuesIn.data(), 0, numValues, deviceStreamIn, transferKind, nullptr);
                copyBetweenDeviceBuffers(
                        &bufferOut, &bufferIn, numValues, deviceStreamIn, transferKind, nullptr);
                if (transferKind == GpuApiCallBehavior::Async)
                {
                    deviceStreamIn.synchronize();
                }
                copyFromDeviceBuffer(
                        valuesOut.data(), &bufferOut, 0, numValues, deviceStreamOut, transferKind, nullptr);
                if (transferKind == GpuApiCallBehavior::Async)
                {
                    deviceStreamOut.synchronize();
                }
                EXPECT_THAT(valuesOut, Pointwise(Eq(), valuesIn))
                        << "Changed after H2D, D2D and D2H " << enumValueToString(transferKind)
                        << " copy.";
                freeDeviceBuffer(&bufferIn);
                freeDeviceBuffer(&bufferOut);
            }
        }
    }
}

#    endif // GMX_GPU_CUDA


} // namespace
} // namespace test
} // namespace gmx

#endif // GMX_GPU
