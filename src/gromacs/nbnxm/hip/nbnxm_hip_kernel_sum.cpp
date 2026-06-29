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
 *  \brief
 *  Explicitly instantiate NBNXM HIP kernel for sum up
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "nbnxm_hip_kernel_sum.h"

#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_hip_kernel_utils.h"
#include "nbnxm_hip_types.h"

namespace gmx
{

namespace
{

__forceinline__ __device__ void float3ReduceFinal(float3* input_ptr)
{
    const unsigned int flat_id = threadIdx.x;
    const unsigned int size    = c_numShiftVectors;
    float3             input;
    input.x = atomicExch(&(input_ptr[size * (flat_id + 1)].x), 0.0F);
    input.y = atomicExch(&(input_ptr[size * (flat_id + 1)].x) + 1, 0.0F);
    input.z = atomicExch(&(input_ptr[size * (flat_id + 1)].x) + 2, 0.0F);

#pragma unroll
    for (unsigned int offset = 1; offset < warpSize; offset *= 2)
    {
        input.x = input.x + __shfl_down(input.x, offset);
        input.y = input.y + __shfl_down(input.y, offset);
        input.z = input.z + __shfl_down(input.z, offset);
    }

    if (flat_id == 0 || flat_id == warpSize)
    {
        atomicAdd(&(input_ptr[0].x), input.x);
        atomicAdd(&(input_ptr[0].x) + 1, input.y);
        atomicAdd(&(input_ptr[0].x) + 2, input.z);
    }
}

__forceinline__ __device__ void energyReduceFinal(float* e_lj_ptr, float* e_el_ptr)
{
    const unsigned int flat_id = threadIdx.x;

    float E_lj = atomicExch(e_lj_ptr + (flat_id + 1), 0.0F);
    float E_el = atomicExch(e_el_ptr + (flat_id + 1), 0.0F);

#pragma unroll
    for (unsigned int offset = 1; offset < warpSize; offset *= 2)
    {
        E_lj += __shfl_down(E_lj, offset);
        E_el += __shfl_down(E_el, offset);
    }

    if (flat_id == 0 || flat_id == warpSize)
    {
        atomicAdd(e_lj_ptr, E_lj);
        atomicAdd(e_el_ptr, E_el);
    }
}

template<unsigned int BlockSize, bool computeEnergy, bool computeVirial>
__launch_bounds__(BlockSize) __global__ void nbnxmKernelSumUp(NBAtomDataGpu atdat)
{
    const unsigned int bidx = blockIdx.x;


    // Sum up fshifts
    if constexpr (computeVirial)
    {
        float3* values_ptr = asFloat3(atdat.fShift) + bidx;
        float3ReduceFinal(values_ptr);
    }

    // Sum up energies
    if constexpr (computeEnergy)
    {
        if (bidx == 0)
        {
            energyReduceFinal(atdat.eLJ, atdat.eElec);
        }
    }
}

//! \brief NBNXM kernel launch code.
template<int blockSize, bool computeEnergy, bool computeVirial, class... Args>
void launchKernelSumUp(const DeviceStream& deviceStream, Args*... args)
{
    KernelLaunchConfig config;
    config.gridSize[0]      = computeVirial ? c_numShiftVectors : 1;
    config.blockSize[0]     = blockSize;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.sharedMemorySize = 0;

    auto kernel = nbnxmKernelSumUp<blockSize, computeEnergy, computeVirial>;

    const auto kernelArgs = prepareGpuKernelArguments(kernel, config, args...);

    launchGpuKernel(kernel, config, deviceStream, nullptr, "nbnxm_hip_sum_up", kernelArgs);
}

template<int blockSize>
void chooseAndLaunchNbnxmKernelSumUp(NbnxmGpu* nb, const StepWorkload& stepWork, const InteractionLocality iloc)
{
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
    NBAtomDataGpu*      adat         = nb->atdat;
    dispatchTemplatedFunction(
            [&](auto computeEnergy_, auto computeVirial_)
            {
                if constexpr (computeEnergy_ || computeVirial_)
                {
                    launchKernelSumUp<blockSize, computeEnergy_, computeVirial_>(deviceStream, adat);
                }
            },
            stepWork.computeEnergy,
            stepWork.computeVirial);
}

} // namespace

void launchNbnxmKernelSumUp(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc)
{
    chooseAndLaunchNbnxmKernelSumUp<sc_energyVirialNumElementsSeparateDeviceReduction>(nb, stepWork, iloc);
}

} // namespace gmx
