/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Runners for tests of CUDA types compatibility.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "typecasts_runner.h"

#include <vector>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace test
{

/* \brief Perform a component-wise conversion of the float3 vector back to RVec format.
 *
 * This is needed to pass the data back to the CPU testing code for comparison with the initial input.
 *
 * \param[out] rVecOutput    Output data in RVec format for the output.
 * \param[in]  float3Output  Output data in float3 format.
 * \param[in]  numElements   Size of the data buffers.
 */
void inline saveFloat3InRVecFormat(ArrayRef<gmx::RVec> rVecOutput, const float3* float3Output, int numElements)
{
    for (int i = 0; i < numElements; i++)
    {
        rVecOutput[i][XX] = float3Output[i].x;
        rVecOutput[i][YY] = float3Output[i].y;
        rVecOutput[i][ZZ] = float3Output[i].z;
    }
}

void convertRVecToFloat3OnHost(ArrayRef<gmx::RVec> rVecOutput, ArrayRef<const gmx::RVec> rVecInput)
{
    const int numElements = rVecInput.size();

    float3* dataFloat3 = asFloat3(const_cast<RVec*>(rVecInput.data()));

    saveFloat3InRVecFormat(rVecOutput, dataFloat3, numElements);
}

//! Number of CUDA threads in a block.
constexpr static int c_threadsPerBlock = 256;

/*! \brief GPU kernel to perform type conversion on the device.
 *
 * \param[out] gm_float3Output Buffer to write the output into.
 * \param[in]  gm_rVecInput    Input data in RVec format.
 * \param[in]  size            Size of the data buffers.
 *
 */
static __global__ void convertRVecToFloat3OnDevice_kernel(DeviceBuffer<float3> gm_float3Output,
                                                          DeviceBuffer<RVec>   gm_rVecInput,
                                                          const int            size)
{
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadIndex < size)
    {
        gm_float3Output[threadIndex] = asFloat3(gm_rVecInput)[threadIndex];
    }
}

void convertRVecToFloat3OnDevice(ArrayRef<gmx::RVec>       h_rVecOutput,
                                 ArrayRef<const gmx::RVec> h_rVecInput,
                                 const TestDevice*         testDevice)
{
    const DeviceContext& deviceContext = testDevice->deviceContext();
    const DeviceStream&  deviceStream  = testDevice->deviceStream();
    deviceContext.activate();

    const int numElements = h_rVecInput.size();

    DeviceBuffer<RVec> d_rVecInput;
    allocateDeviceBuffer(&d_rVecInput, numElements, deviceContext);
    copyToDeviceBuffer(
            &d_rVecInput, h_rVecInput.data(), 0, numElements, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    DeviceBuffer<float3> d_float3Output;
    allocateDeviceBuffer(&d_float3Output, numElements * DIM, deviceContext);

    std::vector<float3> h_float3Output(numElements);

    KernelLaunchConfig kernelLaunchConfig;
    kernelLaunchConfig.gridSize[0]      = gmx::divideRoundUp(numElements, c_threadsPerBlock);
    kernelLaunchConfig.blockSize[0]     = c_threadsPerBlock;
    kernelLaunchConfig.blockSize[1]     = 1;
    kernelLaunchConfig.blockSize[2]     = 1;
    kernelLaunchConfig.sharedMemorySize = 0;

    auto       kernelPtr  = convertRVecToFloat3OnDevice_kernel;
    const auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, kernelLaunchConfig, &d_float3Output, &d_rVecInput, &numElements);
    launchGpuKernel(kernelPtr,
                    kernelLaunchConfig,
                    deviceStream,
                    nullptr,
                    "convertRVecToFloat3OnDevice_kernel",
                    kernelArgs);

    copyFromDeviceBuffer(
            h_float3Output.data(), &d_float3Output, 0, numElements, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    saveFloat3InRVecFormat(h_rVecOutput, h_float3Output.data(), numElements);

    freeDeviceBuffer(&d_rVecInput);
    freeDeviceBuffer(&d_float3Output);
}

} // namespace test
} // namespace gmx
