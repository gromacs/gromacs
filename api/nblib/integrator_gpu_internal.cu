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
 * Implements nblib integrator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include "nblib/integrator_gpu_internal.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/typecasts.cuh"

namespace nblib
{

__global__ void
integrationKernel(real dt, int numParticles, const real* inverseMasses, float4* x, float3* v, float3* f)
{
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < numParticles)
    {
        v[tid].x += f[tid].x * inverseMasses[tid] * dt;
        v[tid].y += f[tid].y * inverseMasses[tid] * dt;
        v[tid].z += f[tid].z * inverseMasses[tid] * dt;

        x[tid].x += v[tid].x * dt;
        x[tid].y += v[tid].y * dt;
        x[tid].z += v[tid].z * dt;
    }
}

void launchLeapFrogKernel(int                  numAtoms,
                          float                dt,
                          DeviceBuffer<Float4> coordinates,
                          DeviceBuffer<Float3> velocities,
                          DeviceBuffer<Float3> forces,
                          DeviceBuffer<float>  inverseMasses,
                          const DeviceStream&  deviceStream)
{
    constexpr unsigned numThreads = 256;
    unsigned           numBlocks  = (numAtoms - 1) / numThreads + 1;

    KernelLaunchConfig launchConfig;
    launchConfig.gridSize[0]  = numBlocks;
    launchConfig.blockSize[0] = numThreads;

    auto       kernelPtr  = integrationKernel;
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      launchConfig,
                                                      &dt,
                                                      &numAtoms,
                                                      &inverseMasses,
                                                      reinterpret_cast<float4**>(&coordinates),
                                                      asFloat3Pointer(&velocities),
                                                      asFloat3Pointer(&forces));
    launchGpuKernel(kernelPtr, launchConfig, deviceStream, nullptr, "leapfrog_kernel", kernelArgs);
}

} // namespace nblib
