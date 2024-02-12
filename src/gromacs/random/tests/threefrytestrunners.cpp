/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief GPU Runners for the ThreeFry random engine
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_random
 */

#include "gmxpre.h"

#include "threefrytestrunners.h"

#include "config.h"

#include <gtest/gtest.h>

#if GPU_THREEFRY_SUPPORTED

#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/random/threefry.h"
#    include "gromacs/utility/exceptions.h"
#    if GMX_GPU_SYCL
#        include "gromacs/gpu_utils/gmxsycl.h"
#    endif

namespace
{

#    if GMX_GPU_CUDA || GMX_GPU_HIP

template<unsigned int rounds, unsigned int internalCounterBits>
GMX_KERNEL_ATTRIBUTE void setupDeviceRngKeys(gmx::ThreeFry2x64General<rounds, internalCounterBits>* d_rng,
                                             uint64_t key0,
                                             uint64_t key1)
{
    new (d_rng) gmx::ThreeFry2x64General<rounds, internalCounterBits>(key0, key1);
}

template<unsigned int rounds, unsigned int internalCounterBits>
GMX_KERNEL_ATTRIBUTE void restart(gmx::ThreeFry2x64General<rounds, internalCounterBits>* d_rng,
                                  uint64_t                                               ctr0,
                                  uint64_t                                               ctr1)
{
    d_rng->restart(ctr0, ctr1);
}

template<unsigned int rounds, unsigned int internalCounterBits>
GMX_KERNEL_ATTRIBUTE void next(gmx::ThreeFry2x64General<rounds, internalCounterBits>* d_rng,
                               uint64_t* __restrict__ gm_result)
{
    *gm_result = (*d_rng)();
}

#    elif GMX_GPU_SYCL

// SYCL kernel-name tags for the three ThreeFry test kernels.
template<unsigned int rounds, unsigned int internalCounterBits>
class ThreeFrySetupKernel;
template<unsigned int rounds, unsigned int internalCounterBits>
class ThreeFryRestartKernel;
template<unsigned int rounds, unsigned int internalCounterBits>
class ThreeFryNextKernel;

#    endif

} // namespace

namespace gmx
{
namespace test
{

template<unsigned int rounds, unsigned int internalCounterBits>
ThreeFryGeneralDeviceTestRunner<rounds, internalCounterBits>::ThreeFryGeneralDeviceTestRunner(
        const TestDevice& testDevice,
        uint64_t          key0,
        uint64_t          key1) :
    testDevice_(testDevice)
{
    const DeviceStream&  deviceStream  = testDevice_.deviceStream();
    const DeviceContext& deviceContext = testDevice_.deviceContext();
    deviceContext.activate();
    allocateDeviceBuffer(&rng_, 1, deviceContext);
    allocateDeviceBuffer(&d_result, 1, deviceContext);

#    if GMX_GPU_CUDA || GMX_GPU_HIP
    KernelLaunchConfig kernelLaunchConfig;

    auto kernelPtr = setupDeviceRngKeys<rounds, internalCounterBits>;
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &rng_, &key0, &key1);

    launchGpuKernel(kernelPtr,
                    kernelLaunchConfig,
                    deviceStream,
                    nullptr,
                    "test_threefry_setup_device_rng_keys_kernel",
                    kernelArgs);
#    elif GMX_GPU_SYCL
    using Rng   = gmx::ThreeFry2x64General<rounds, internalCounterBits>;
    Rng* gm_rng = rng_.get_pointer();
    deviceStream.stream().submit(
            [&](sycl::handler& cgh)
            {
                cgh.single_task<ThreeFrySetupKernel<rounds, internalCounterBits>>(
                        [=]() { new (gm_rng) Rng(key0, key1); });
            });
#    endif
}

template<unsigned int rounds, unsigned int internalCounterBits>
ThreeFryGeneralDeviceTestRunner<rounds, internalCounterBits>::~ThreeFryGeneralDeviceTestRunner()
{
    freeDeviceBuffer(&rng_);
    freeDeviceBuffer(&d_result);
}

template<unsigned int rounds, unsigned int internalCounterBits>
void ThreeFryGeneralDeviceTestRunner<rounds, internalCounterBits>::restartRng(uint64_t ctr0, uint64_t ctr1)
{
    const DeviceStream&  deviceStream  = testDevice_.deviceStream();
    const DeviceContext& deviceContext = testDevice_.deviceContext();
    deviceContext.activate();

#    if GMX_GPU_CUDA || GMX_GPU_HIP
    KernelLaunchConfig kernelLaunchConfig;

    auto kernelPtr = restart<rounds, internalCounterBits>;
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &rng_, &ctr0, &ctr1);

    launchGpuKernel(
            kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "test_threefry_restart_kernel", kernelArgs);
#    elif GMX_GPU_SYCL
    using Rng   = gmx::ThreeFry2x64General<rounds, internalCounterBits>;
    Rng* gm_rng = rng_.get_pointer();
    deviceStream.stream().submit(
            [&](sycl::handler& cgh)
            {
                cgh.single_task<ThreeFryRestartKernel<rounds, internalCounterBits>>(
                        [=]() { gm_rng->restart(ctr0, ctr1); });
            });
#    endif
}

template<unsigned int rounds, unsigned int internalCounterBits>
uint64_t ThreeFryGeneralDeviceTestRunner<rounds, internalCounterBits>::nextRng()
{
    const DeviceStream& deviceStream = testDevice_.deviceStream();

#    if GMX_GPU_CUDA || GMX_GPU_HIP
    KernelLaunchConfig kernelLaunchConfig;

    auto kernelPtr = next<rounds, internalCounterBits>;
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &rng_, &d_result);
    launchGpuKernel(
            kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "test_threefry_next_kernel", kernelArgs);
#    elif GMX_GPU_SYCL
    using Rng         = gmx::ThreeFry2x64General<rounds, internalCounterBits>;
    Rng*      gm_rng  = rng_.get_pointer();
    uint64_t* gm_dest = d_result.get_pointer();
    deviceStream.stream().submit(
            [&](sycl::handler& cgh)
            {
                cgh.single_task<ThreeFryNextKernel<rounds, internalCounterBits>>(
                        [=]() { *gm_dest = (*gm_rng)(); });
            });
#    endif

    uint64_t h_result;
    copyFromDeviceBuffer(&h_result, &d_result, 0, 1, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    return h_result;
}

} // namespace test
} // namespace gmx

template class gmx::test::ThreeFryGeneralDeviceTestRunner<13, 0>;
template class gmx::test::ThreeFryGeneralDeviceTestRunner<20, 0>;
template class gmx::test::ThreeFryGeneralDeviceTestRunner<40, 0>;

#endif
