/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * NVSHMEM Put/Allocator/Init tests of NVSHMEM support.
 *
 * \author Mahesh Doijade <mdoijade@nvidia.com>
 */
#include "gmxpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/gmxmpi.h"

#if GMX_NVSHMEM
#    include <nvshmem.h>
#    include <nvshmemx.h>

#    include "gromacs/gpu_utils/devicebuffer.cuh"
#    include "gromacs/hardware/device_information.h"
#    include "gromacs/utility/exceptions.h"
#    include "gromacs/utility/stringutil.h"

#    include "testutils/test_device.h"

#    include "nvshmem_simple_put.cuh"

namespace gmx
{

namespace test
{

//! Number of CUDA threads in a block.
constexpr static int c_threadsPerBlock = 1;

/*! \brief GPU kernel to perform nvshmem puts on symmetric buffer.
 *
 * \param[out] gm_outDestination Buffer to write the output into.
 * \param[in]  size            Size of the data buffers.
 *
 */
static __global__ void nvshmemSimplePut(DeviceBuffer<int> gm_outDestination, const int size)
{
// NVSHMEM is supported on Volta+
#    if __CUDA_ARCH__ >= 700
    int mype = nvshmem_my_pe();
    int npes = nvshmem_n_pes();
    int peer = (mype + 1) % npes;

    if (threadIdx.x < size)
    {
        nvshmem_int_p(gm_outDestination + threadIdx.x, mype, peer);
    }
#    endif
}

void nvshmemRunSimplePut(const TestDevice* testDevice)
{
    const DeviceContext& deviceContext = testDevice->deviceContext();
    const DeviceStream&  deviceStream  = testDevice->deviceStream();
    deviceContext.activate();

    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    nvshmemx_init_attr_t attr;
    attr.mpi_comm = &mpi_comm;
    nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
    int mype_node = nvshmem_team_my_pe(NVSHMEMX_TEAM_NODE);

    const int         numElements = 1;
    DeviceBuffer<int> destination;
    allocateDeviceBufferNvShmem(&destination, 1, deviceContext);

    KernelLaunchConfig kernelLaunchConfig;
    kernelLaunchConfig.gridSize[0]      = gmx::divideRoundUp(numElements, c_threadsPerBlock);
    kernelLaunchConfig.blockSize[0]     = c_threadsPerBlock;
    kernelLaunchConfig.blockSize[1]     = 1;
    kernelLaunchConfig.blockSize[2]     = 1;
    kernelLaunchConfig.sharedMemorySize = 0;

    auto       kernelPtr = nvshmemSimplePut;
    const auto kernelArgs =
            prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig, &destination, &numElements);
    launchGpuKernel(kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "nvshmemSimplePut", kernelArgs);

    nvshmemx_barrier_all_on_stream(deviceStream.stream());
    int msg = -1;

    copyFromDeviceBuffer(&msg, &destination, 0, numElements, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    freeDeviceBuffer(&destination);
    nvshmem_finalize();
    EXPECT_EQ(msg, (mype_node + 1) % 2);
}

} // namespace test
} // namespace gmx

#endif
