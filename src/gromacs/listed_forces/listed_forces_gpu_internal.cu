/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 *
 * \brief Implements CUDA bonded functionality
 *
 * \author Jon Vincent <jvincent@nvidia.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */

#include "gmxpre.h"

#include <math_constants.h>

#include <cassert>
#include <cstddef>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops_cuda.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/listed_forces/listed_forces_gpu_internal_shared.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/gmxassert.h"

#include "listed_forces_gpu_impl.h"

#if defined(_MSVC)
#    include <limits>
#endif


namespace gmx
{

template<bool calcVir, bool calcEner>
__global__ void bonded_kernel_gpu(BondedGpuKernelParameters kernelParams,
                                  BondedGpuKernelBuffers    kernelBuffers,
                                  float4*                   gm_xq_in,
                                  float3*                   gm_f,
                                  float3*                   gm_fShift)
{
    GMX_DEVICE_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    const int tid          = blockIdx.x * blockDim.x + threadIdx.x;
    float     vtot_loc     = 0.0F;
    float     vtotElec_loc = 0.0F; // Used only for InteractionFunction::LennardJones14

    extern __shared__ char sm_dynamicShmem[];
    char*                  sm_nextSlotPtr = sm_dynamicShmem;
    float3*                sm_fShiftLoc   = reinterpret_cast<float3*>(sm_nextSlotPtr);
    sm_nextSlotPtr += c_numShiftVectors * sizeof(float3);

    static_assert(sizeof(DeviceFloat4) == sizeof(float4));
    DeviceFloat4* gm_xq = reinterpret_cast<DeviceFloat4*>(gm_xq_in);

    if (calcVir)
    {
        if (threadIdx.x < c_numShiftVectors)
        {
            sm_fShiftLoc[threadIdx.x] = make_float3(0.0F, 0.0F, 0.0F);
        }
        __syncthreads();
    }

    InteractionFunction fType;
    bool                threadComputedPotential = false;
#pragma unroll
    for (int j = 0; j < numFTypesOnGpu; j++)
    {
        if (tid >= kernelParams.fTypeRangeStart[j] && tid <= kernelParams.fTypeRangeEnd[j])
        {
            const int      numBonds = kernelParams.numFTypeBonds[j];
            int            fTypeTid = tid - kernelParams.fTypeRangeStart[j];
            const t_iatom* iatoms   = kernelBuffers.d_iatoms[j];
            fType                   = kernelParams.fTypesOnGpu[j];
            if (calcEner)
            {
                threadComputedPotential = true;
            }

            if (fTypeTid >= numBonds)
            {
                break;
            }


            switch (fType)
            {
                case InteractionFunction::Bonds:
                    bonds_gpu<calcVir, calcEner>(fTypeTid,
                                                 &vtot_loc,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc,
                                                 threadIdx.x);
                    break;
                case InteractionFunction::Angles:
                    angles_gpu<calcVir, calcEner>(fTypeTid,
                                                  &vtot_loc,
                                                  iatoms,
                                                  kernelBuffers.d_forceParams,
                                                  gm_xq,
                                                  gm_f,
                                                  sm_fShiftLoc,
                                                  kernelParams.pbcAiuc,
                                                  threadIdx.x);
                    break;
                case InteractionFunction::UreyBradleyPotential:
                    urey_bradley_gpu<calcVir, calcEner>(fTypeTid,
                                                        &vtot_loc,
                                                        iatoms,
                                                        kernelBuffers.d_forceParams,
                                                        gm_xq,
                                                        gm_f,
                                                        sm_fShiftLoc,
                                                        kernelParams.pbcAiuc,
                                                        threadIdx.x);
                    break;
                case InteractionFunction::ProperDihedrals:
                case InteractionFunction::PeriodicImproperDihedrals:
                    pdihs_gpu<calcVir, calcEner>(fTypeTid,
                                                 &vtot_loc,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc,
                                                 threadIdx.x);
                    break;
                case InteractionFunction::RyckaertBellemansDihedrals:
                    rbdihs_gpu<calcVir, calcEner>(fTypeTid,
                                                  &vtot_loc,
                                                  iatoms,
                                                  kernelBuffers.d_forceParams,
                                                  gm_xq,
                                                  gm_f,
                                                  sm_fShiftLoc,
                                                  kernelParams.pbcAiuc,
                                                  threadIdx.x);
                    break;
                case InteractionFunction::ImproperDihedrals:
                    idihs_gpu<calcVir, calcEner>(fTypeTid,
                                                 &vtot_loc,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc,
                                                 threadIdx.x);
                    break;
                case InteractionFunction::LennardJones14:
                    pairs_gpu<calcVir, calcEner>(fTypeTid,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc,
                                                 kernelParams.electrostaticsScaleFactor,
                                                 &vtot_loc,
                                                 &vtotElec_loc,
                                                 threadIdx.x);
                    break;
                default:
                    // these types do not appear on the GPU
                    break;
            }
            break;
        }
    }

    if (threadComputedPotential)
    {
        float* vtot = kernelBuffers.d_vTot + static_cast<ptrdiff_t>(fType);
        float* vtotElec = kernelBuffers.d_vTot + static_cast<ptrdiff_t>(InteractionFunction::Coulomb14);

        // Perform warp-local reduction
        vtot_loc += __shfl_down_sync(c_fullWarpMask, vtot_loc, 1);
        vtotElec_loc += __shfl_up_sync(c_fullWarpMask, vtotElec_loc, 1);
        if (threadIdx.x & 1)
        {
            vtot_loc = vtotElec_loc;
        }
#pragma unroll 4
        for (int i = 2; i < warpSize; i *= 2)
        {
            vtot_loc += __shfl_down_sync(c_fullWarpMask, vtot_loc, i);
        }

        // Write reduced warp results into global memory
        if (threadIdx.x % warpSize == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        else if ((threadIdx.x % warpSize == 1) && (fType == InteractionFunction::LennardJones14))
        {
            atomicAdd(vtotElec, vtot_loc);
        }
    }
    /* Accumulate shift vectors from shared memory to global memory on the first c_numShiftVectors threads of the block. */
    if (calcVir)
    {
        __syncthreads();
        if (threadIdx.x < c_numShiftVectors)
        {
            staggeredAtomicAddForce(&gm_fShift[threadIdx.x], sm_fShiftLoc[threadIdx.x], threadIdx.x);
        }
    }
}


/*-------------------------------- End CUDA kernels-----------------------------*/


template<bool calcVir, bool calcEner>
void ListedForcesGpu::Impl::launchKernel()
{
    GMX_ASSERT(haveInteractions_,
               "Cannot launch bonded GPU kernels unless bonded GPU work was scheduled");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuBonded);

    int fTypeRangeEnd = kernelParams_.fTypeRangeEnd[numFTypesOnGpu - 1];

    if (fTypeRangeEnd < 0)
    {
        return;
    }

    auto kernelPtr = bonded_kernel_gpu<calcVir, calcEner>;

    const auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, kernelLaunchConfig_, &kernelParams_, &kernelBuffers_, &d_xq_, &d_f_, &d_fShift_);

    launchGpuKernel(kernelPtr,
                    kernelLaunchConfig_,
                    deviceStream_,
                    nullptr,
                    "bonded_kernel_gpu<calcVir, calcEner>",
                    kernelArgs);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuBonded);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void ListedForcesGpu::launchKernel(const gmx::StepWorkload& stepWork)
{
    if (stepWork.computeEnergy)
    {
        // When we need the energy, we also need the virial
        impl_->launchKernel<true, true>();
    }
    else if (stepWork.computeVirial)
    {
        impl_->launchKernel<true, false>();
    }
    else
    {
        impl_->launchKernel<false, false>();
    }
}

} // namespace gmx
