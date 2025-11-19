/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * \brief Implements SYCL GPU bonded functionality
 *
 * \author Andrey Alekseenko <al42and@gmail.com>
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

#include "config.h"

#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/listed_forces/listed_forces_gpu_internal_shared.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc_aiuc_sycl.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"

#include "listed_forces_gpu_impl.h"

template<bool calcVir, bool calcEner>
class BondedKernel;

namespace gmx
{

using sycl::access::fence_space;
using mode = sycl::access_mode;

template<bool calcVir, bool calcEner, typename CommandGroupHandler>
auto bondedKernel(CommandGroupHandler              cgh,
                  const BondedGpuKernelParameters& kernelParams,
                  const DeviceBuffer<t_iatom>      gm_iatoms_[numFTypesOnGpu],
                  float* __restrict__ gm_vTot,
                  const t_iparams* __restrict__ gm_forceParams_,
                  const sycl::float4* __restrict__ gm_xq_in_,
                  Float3* __restrict__ gm_f_,
                  Float3* __restrict__ gm_fShift_)
{
    static_assert(sizeof(DeviceFloat4) == sizeof(sycl::float4));
    const DeviceFloat4*             gm_xq_ = reinterpret_cast<const DeviceFloat4*>(gm_xq_in_);
    sycl::global_ptr<const t_iatom> gm_iatomsTemp[numFTypesOnGpu];
    for (int i = 0; i < numFTypesOnGpu; i++)
    {
        gm_iatomsTemp[i] = gm_iatoms_[i].get_pointer();
    }
    const FTypeArray<sycl::global_ptr<const t_iatom>> gm_iatoms(gm_iatomsTemp);

    const FTypeArray<int> fTypeRangeStart(kernelParams.fTypeRangeStart);
    const FTypeArray<int> fTypeRangeEnd(kernelParams.fTypeRangeEnd);
    const FTypeArray<int> numFTypeBonds(kernelParams.numFTypeBonds);

    const auto electrostaticsScaleFactor = kernelParams.electrostaticsScaleFactor;

    using FShiftLoc              = StaticLocalStorage<Float3, c_numShiftVectors>;
    auto sm_fShiftLocHostStorage = FShiftLoc::makeHostStorage(cgh);

    const PbcAiuc pbcAiuc = kernelParams.pbcAiuc;

    return [=](sycl::nd_item<1> itemIdx)
    {
        // This declaration works on the device
        typename FShiftLoc::DeviceStorage sm_fShiftLocDeviceStorage;
        // Extract the valid pointer to local storage
        sycl::local_ptr<Float3> sm_fShiftLoc =
                FShiftLoc::get_pointer(sm_fShiftLocHostStorage, sm_fShiftLocDeviceStorage);

        sycl::global_ptr<const t_iparams>    gm_forceParams = gm_forceParams_;
        sycl::global_ptr<const DeviceFloat4> gm_xq          = gm_xq_;
        sycl::global_ptr<Float3>             gm_f           = gm_f_;
        sycl::global_ptr<Float3>             gm_fShift      = gm_fShift_;

        const int tid          = itemIdx.get_global_linear_id();
        const int localId      = itemIdx.get_local_linear_id();
        float     vtot_loc     = 0.0F;
        float     vtotElec_loc = 0.0F;

        if constexpr (calcVir)
        {
            if (localId < c_numShiftVectors)
            {
                sm_fShiftLoc[localId] = { 0.0F, 0.0F, 0.0F };
            }
            itemIdx.barrier(fence_space::local_space);
        }

        InteractionFunction fType;
        bool                threadComputedPotential = false;
#pragma unroll
        for (int j = 0; j < numFTypesOnGpu; j++)
        {
            if (tid >= fTypeRangeStart[j] && tid <= fTypeRangeEnd[j])
            {
                const int                             numBonds = numFTypeBonds[j];
                const int                             fTypeTid = tid - fTypeRangeStart[j];
                const sycl::global_ptr<const t_iatom> iatoms   = gm_iatoms[j];
                fType                                          = fTypesOnGpu[j];

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
                        bonds_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case InteractionFunction::Angles:
                        angles_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case InteractionFunction::UreyBradleyPotential:
                        urey_bradley_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case InteractionFunction::ProperDihedrals:
                    case InteractionFunction::PeriodicImproperDihedrals:
                        pdihs_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case InteractionFunction::RyckaertBellemansDihedrals:
                        rbdihs_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case InteractionFunction::ImproperDihedrals:
                        idihs_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case InteractionFunction::LennardJones14:
                        pairs_gpu<calcVir, calcEner>(fTypeTid,
                                                     iatoms,
                                                     gm_forceParams,
                                                     gm_xq,
                                                     gm_f,
                                                     sm_fShiftLoc,
                                                     pbcAiuc,
                                                     electrostaticsScaleFactor,
                                                     &vtot_loc,
                                                     &vtotElec_loc,
                                                     localId);
                        break;
                    default:
                        // these types do not appear on the GPU
                        break;
                }
                break;
            }
        }

        if (calcEner && threadComputedPotential)
        {
            subGroupBarrier(itemIdx); // Should not be needed, but https://github.com/AdaptiveCpp/AdaptiveCpp/issues/823
            sycl::sub_group sg = itemIdx.get_sub_group();
            vtot_loc           = sycl::reduce_over_group(sg, vtot_loc, sycl::plus<float>());
            vtotElec_loc       = sycl::reduce_over_group(sg, vtotElec_loc, sycl::plus<float>());
            if (sg.leader())
            {
                atomicFetchAdd(gm_vTot[static_cast<int>(fType)], vtot_loc);
                if (fType == InteractionFunction::LennardJones14)
                {
                    atomicFetchAdd(gm_vTot[static_cast<int>(InteractionFunction::Coulomb14)], vtotElec_loc);
                }
            }
        }
        /* Accumulate shift vectors from shared memory to global memory on the first c_numShiftVectors threads of the block. */
        if constexpr (calcVir)
        {
            itemIdx.barrier(fence_space::local_space);
            if (localId < c_numShiftVectors)
            {
                const Float3 tmp = sm_fShiftLoc[localId];
                atomicFetchAdd(gm_fShift[localId], tmp);
            }
        }
    };
}


template<bool calcVir, bool calcEner>
void ListedForcesGpu::Impl::launchKernel()
{
    GMX_ASSERT(haveInteractions_,
               "Cannot launch bonded GPU kernels unless bonded GPU work was scheduled");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuBonded);

    if (kernelParams_.fTypeRangeEnd[numFTypesOnGpu - 1] < 0)
    {
        return;
    }

    using kernelNameType = BondedKernel<calcVir, calcEner>;

    const sycl::nd_range<1> rangeAll(kernelLaunchConfig_.blockSize[0] * kernelLaunchConfig_.gridSize[0],
                                     kernelLaunchConfig_.blockSize[0]);

    auto kernelFunctionBuilder = bondedKernel<calcVir, calcEner, CommandGroupHandler>;
    syclSubmitWithoutEvent<kernelNameType>(deviceStream_.stream(),
                                           kernelFunctionBuilder,
                                           rangeAll,
                                           kernelParams_,
                                           kernelBuffers_.d_iatoms,
                                           kernelBuffers_.d_vTot.get_pointer(),
                                           kernelBuffers_.d_forceParams.get_pointer(),
                                           d_xq_.get_pointer(),
                                           d_f_.get_pointer(),
                                           d_fShift_.get_pointer());

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
