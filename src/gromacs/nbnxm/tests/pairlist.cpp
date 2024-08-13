/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for gpu pairlist datastructure memory handling.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "config.h"

#include "gmock/gmock.h"

#if GMX_GPU && !GMX_GPU_HIP

#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gpu_utils.h"
#    include "gromacs/gpu_utils/hostallocator.h"
#    include "gromacs/hardware/device_management.h"
#    include "gromacs/nbnxm/gpu_types_common.h"

#    include "testutils/test_hardware_environment.h"
#    include "testutils/testasserts.h"

using ::testing::Eq;
using ::testing::Pointwise;


namespace gmx
{

namespace test
{

namespace
{

//! Helper wrapper for host resident data
struct HostBuffers
{
    HostVector<nbnxn_sci_t>       h_sci;
    HostVector<nbnxn_cj_packed_t> h_cjPacked;
    HostVector<nbnxn_excl_t>      h_excl;
    HostVector<int>               h_rollingPrunePart;

    //! We only build this with our test data size and with pinning when possible
    HostBuffers(int dataSize)
    {
        changePinningPolicy(&h_sci, gmx::PinningPolicy::PinnedIfSupported);
        changePinningPolicy(&h_cjPacked, gmx::PinningPolicy::PinnedIfSupported);
        changePinningPolicy(&h_excl, gmx::PinningPolicy::PinnedIfSupported);
        changePinningPolicy(&h_rollingPrunePart, gmx::PinningPolicy::PinnedIfSupported);

        h_sci.resize(dataSize);
        h_cjPacked.resize(dataSize);
        h_excl.resize(dataSize);
        h_rollingPrunePart.resize(dataSize);
    }
};

//! Allocate device side buffers for testing
void allocateDeviceBuffers(GpuPairlist* pairlist, int allocationSize, const DeviceContext& deviceContext)
{
    pairlist->sciAllocationSize                  = allocationSize;
    pairlist->packedJClustersAllocationSize      = allocationSize;
    pairlist->exclAllocationSize                 = allocationSize;
    pairlist->d_rollingPruningPartAllocationSize = allocationSize;

    allocateDeviceBuffer(&pairlist->sci, allocationSize, deviceContext);
    allocateDeviceBuffer(&pairlist->cjPacked, allocationSize, deviceContext);
    allocateDeviceBuffer(&pairlist->excl, allocationSize, deviceContext);
    allocateDeviceBuffer(&pairlist->d_rollingPruningPart, allocationSize, deviceContext);
}

//! Set up host resident data structures with junk test data.
HostBuffers prepareHostBuffers(int dataSize)
{
    HostBuffers hostbuffers(dataSize);

    int counter = 0;
    for (auto& value : hostbuffers.h_sci)
    {
        value.sci           = counter++;
        value.shift         = counter++;
        value.cjPackedBegin = counter++;
        value.cjPackedEnd   = counter++;
    }

    for (auto& value : hostbuffers.h_cjPacked)
    {
        for (int index = 0; index < c_nbnxnGpuClusterpairSplit; index++)
        {
            value.imei[index].excl_ind = counter++;
        }
        for (int index = 0; index < c_nbnxnGpuJgroupSize; index++)
        {
            value.cj[index] = counter++;
        }
    }

    for (auto& value : hostbuffers.h_excl)
    {
        for (int index = 0; index < c_nbnxnGpuExclSize; index++)
        {
            value.pair[index] = counter++;
        }
    }

    for (auto& value : hostbuffers.h_rollingPrunePart)
    {
        value = counter++;
    }

    return hostbuffers;
}

//! H2D data transfer. Only for checking device side data, transfer is tested elsewhere
void transferHostToDevice(GpuPairlist*        pairlist,
                          const HostBuffers&  hostbuffers,
                          int                 dataSize,
                          const DeviceStream& deviceStream)
{

    copyToDeviceBuffer(
            &pairlist->sci, hostbuffers.h_sci.data(), 0, dataSize, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&pairlist->cjPacked,
                       hostbuffers.h_cjPacked.data(),
                       0,
                       dataSize,
                       deviceStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    copyToDeviceBuffer(
            &pairlist->excl, hostbuffers.h_excl.data(), 0, dataSize, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&pairlist->d_rollingPruningPart,
                       hostbuffers.h_rollingPrunePart.data(),
                       0,
                       dataSize,
                       deviceStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);
}

//! D2H data transfer. Only for checking device side data, transfer is tested elsewhere
void transferDeviceToHost(GpuPairlist* pairlist, HostBuffers* hostbuffers, int dataSize, const DeviceStream& deviceStream)
{

    copyFromDeviceBuffer(
            hostbuffers->h_sci.data(), &pairlist->sci, 0, dataSize, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(hostbuffers->h_cjPacked.data(),
                         &pairlist->cjPacked,
                         0,
                         dataSize,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);
    copyFromDeviceBuffer(
            hostbuffers->h_excl.data(), &pairlist->excl, 0, dataSize, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(hostbuffers->h_rollingPrunePart.data(),
                         &pairlist->d_rollingPruningPart,
                         0,
                         dataSize,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);
}

//! Set up pairlist data structure and transfer data to it
void initPairlistWithData(int                  dataSize,
                          int                  allocationSize,
                          GpuPairlist*         pairlist,
                          const HostBuffers&   hostbuffers,
                          const DeviceContext& deviceContext,
                          const DeviceStream&  deviceStream)
{
    allocateDeviceBuffers(pairlist, allocationSize, deviceContext);

    transferHostToDevice(pairlist, hostbuffers, dataSize, deviceStream);
}

//! Check that all fields in pairlist are properly allocated
void checkPairlistConsistency(const GpuPairlist& pairlist, int allocationSize)
{
    checkDeviceBuffer(pairlist.sci, allocationSize);
    checkDeviceBuffer(pairlist.cjPacked, allocationSize);
    checkDeviceBuffer(pairlist.excl, allocationSize);
    checkDeviceBuffer(pairlist.d_rollingPruningPart, allocationSize);
}

//! Check that device buffer data is valid after transfer to back and forth
void checkValuesFromPairlist(GpuPairlist*        pairlist,
                             const HostBuffers&  inputBuffers,
                             int                 dataSize,
                             const DeviceStream& deviceStream)
{
    HostBuffers fromDevice(dataSize);

    transferDeviceToHost(pairlist, &fromDevice, dataSize, deviceStream);

    EXPECT_THAT(inputBuffers.h_sci, Pointwise(Eq(), fromDevice.h_sci));
    EXPECT_THAT(inputBuffers.h_cjPacked, Pointwise(Eq(), fromDevice.h_cjPacked));
    EXPECT_THAT(inputBuffers.h_excl, Pointwise(Eq(), fromDevice.h_excl));
    EXPECT_THAT(inputBuffers.h_rollingPrunePart, Pointwise(Eq(), fromDevice.h_rollingPrunePart));
}


TEST(GpuPairlistTest, PairlistInitWorks)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        testDevice->deviceContext().activate();
        constexpr int dataSize       = 23;
        constexpr int allocationSize = 42;
        GpuPairlist   pairlist{};
        auto          hostbuffers = prepareHostBuffers(dataSize);
        initPairlistWithData(dataSize,
                             allocationSize,
                             &pairlist,
                             hostbuffers,
                             testDevice->deviceContext(),
                             testDevice->deviceStream());
        checkPairlistConsistency(pairlist, allocationSize);
        checkValuesFromPairlist(&pairlist, hostbuffers, dataSize, testDevice->deviceStream());
    }
}

} // namespace
} // namespace test
} // namespace gmx

#endif
