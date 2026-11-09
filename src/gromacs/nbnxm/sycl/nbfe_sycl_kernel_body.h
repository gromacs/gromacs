/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 *  SYCL nonbonded free-energy kernels for current and foreign lambdas.
 *
 *  This is the SYCL counterpart of the CUDA kernel in
 *  \c nbnxm/cuda/nbfe_cuda_kernel.cuh. As in the CUDA implementation, one
 *  sub-group (warp) processes one entry of the FEP i-atom list and each
 *  work-item of that sub-group processes one j-atom at a time.
 *
 *  \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_BODY_H
#define GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_BODY_H

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/functions.h"
#include "gromacs/utility/template_mp.h"

#include "nbfe_sycl_kernel.h"
#include "nbfe_sycl_kernel_utils.h"
#include "nbnxm_sycl_types.h"

namespace gmx
{

//! \brief Class name for the NBNXM free-energy kernel
template<bool computeForeignLambda, bool doCalcEnergies, bool useSoftCore, enum ElecType elecType, enum VdwType vdwType, int subGroupSize>
class NbfeKernel;

/*! \brief Current-lambda forces/energies or foreign-lambda energies. */
template<int subGroupSize, bool computeForeignLambda, bool doCalcEnergies, bool useSoftCore, enum ElecType elecType, enum VdwType vdwType>
static auto nbfeKernel(const Float4* __restrict__ gm_xq,
                       Float3* __restrict__ gm_f,
                       const Float3* __restrict__ gm_shiftVec,
                       Float3* __restrict__ gm_fShift,
                       float* __restrict__ gm_energyElec,
                       float* __restrict__ gm_energyVdw,
                       float* __restrict__ gm_dvdlElec,
                       float* __restrict__ gm_dvdlVdw,
                       const Float4* __restrict__ gm_q4,
                       const Int4* __restrict__ gm_atomTypes4 /* used iff !ljComb<vdwType> */,
                       const Float4* __restrict__ gm_ljComb4 /* used iff ljComb<vdwType> */,
                       const Float2* __restrict__ gm_nbfp /* used iff !ljComb<vdwType> */,
                       const float* __restrict__ gm_coulombTab /* used iff elecEwaldTab<elecType> */,
                       const int* __restrict__ gm_iinr,
                       const int* __restrict__ gm_jIndex,
                       const int* __restrict__ gm_jjnr,
                       const int* __restrict__ gm_shift,
                       const int* __restrict__ gm_exclFep,
                       const int   numiAtoms,
                       const int   numTypes,
                       const float epsFac,
                       const int   lambdaPower,
                       const float lambdaCoul,
                       const float lambdaVdw,
                       const float* __restrict__ gm_allLambdaCoul,
                       const float* __restrict__ gm_allLambdaVdw,
                       const int                   nLambda,
                       const NbfeInteractionParams params,
                       const bool                  calcShift)
{
    static constexpr EnergyFunctionProperties<elecType, vdwType> props;

    return [=](sycl::nd_item<1> itemIdx) [[sycl::reqd_sub_group_size(subGroupSize)]]
    {
        if constexpr (skipKernelCompilation<subGroupSize>())
        {
            return;
        }

        const sycl::sub_group sg = itemIdx.get_sub_group();
        // Note that the i-list entry index is sub-group-uniform, so the early return
        // below and the sub-group collectives at the end of the kernel are compatible.
        const int tidx    = itemIdx.get_local_linear_id();
        const int laneIdx = tidx % subGroupSize;
        const int iEntry = itemIdx.get_group_linear_id() * sc_nbfeSubGroupsPerGroup + tidx / subGroupSize;

        if (iEntry >= numiAtoms)
        {
            return;
        }

        const float rCutoffMaxSq = props.elecEwaldTwin ? sycl::max(params.rCoulombSq, params.rVdwSq)
                                                       : params.rCoulombSq;

        const int shiftIdx = gm_shift[iEntry];
        const int nj0      = gm_jIndex[iEntry];
        const int nj1      = gm_jIndex[iEntry + 1];
        const int ai       = gm_iinr[iEntry];

        const Float4 xqi = gm_xq[ai];
        const Float3 xi  = Float3(xqi[0], xqi[1], xqi[2]) + gm_shiftVec[shiftIdx];

        const Float4 qi  = gm_q4[ai];
        const float  qAi = qi[0] * epsFac;
        const float  qBi = qi[1] * epsFac;

        int    typeiAB[2] = { 0, 0 };
        Float2 ljcpiAB[2] = { { 0.0F, 0.0F }, { 0.0F, 0.0F } };
        if constexpr (!props.vdwComb)
        {
            const Int4 types = gm_atomTypes4[ai];
            typeiAB[0]       = types[0];
            typeiAB[1]       = types[1];
        }
        else
        {
            const Float4 ljc = gm_ljComb4[ai];
            ljcpiAB[0]       = Float2(ljc[0], ljc[1]);
            ljcpiAB[1]       = Float2(ljc[2], ljc[3]);
        }

        // A uniform lambda loop keeps every lane participating in each reduction,
        // including lanes without a j-atom in the last block.
        for (int lambdaIdx = 0; lambdaIdx <= (computeForeignLambda ? nLambda : 0); lambdaIdx++)
        {
            const float lambdaCoulThis = lambdaIdx == 0 ? lambdaCoul : gm_allLambdaCoul[lambdaIdx - 1];
            const float lambdaVdwThis = lambdaIdx == 0 ? lambdaVdw : gm_allLambdaVdw[lambdaIdx - 1];
            const auto  lf = makeNbfeLambdaFactors(lambdaCoulThis, lambdaVdwThis, lambdaPower);

            Float3           fci = { 0.0F, 0.0F, 0.0F };
            NbfeAccumulators totals;

            for (int i = nj0; i < nj1; i += subGroupSize)
            {
                const int j = i + laneIdx;
                if (j < nj1)
                {
                    const int aj = gm_jjnr[j];

                    NbfePairData pair;
                    pair.pairIncluded = (gm_exclFep == nullptr || gm_exclFep[j] != 0);
                    pair.selfPair     = (ai == aj);
                    // Only updated by nbfeComputeLjParams() when the soft-core is in use
                    pair.alphaCoulombEff = params.alphaCoulomb;
                    pair.alphaVdwEff     = params.alphaVdw;
                    pair.sigma6[0] = pair.sigma6[1] = 0.0F;

                    /* Load j-atom data */
                    const Float4 xqj = gm_xq[aj];
                    const Float3 xj  = Float3(xqj[0], xqj[1], xqj[2]);
                    const Float4 qj  = gm_q4[aj];
                    pair.qq[0]       = qAi * qj[0];
                    pair.qq[1]       = qBi * qj[1];

                    int    typejAB[2] = { 0, 0 };
                    Float2 ljcpjAB[2] = { { 0.0F, 0.0F }, { 0.0F, 0.0F } };
                    if constexpr (!props.vdwComb)
                    {
                        const Int4 types = gm_atomTypes4[aj];
                        typejAB[0]       = types[0];
                        typejAB[1]       = types[1];
                    }
                    else
                    {
                        const Float4 ljc = gm_ljComb4[aj];
                        ljcpjAB[0]       = Float2(ljc[0], ljc[1]);
                        ljcpjAB[1]       = Float2(ljc[2], ljc[3]);
                    }

                    /* Distance between the i- and j-atom */
                    const Float3 rv          = xi - xj;
                    const float  r2Unclamped = rv.norm2();
                    pair.withinCutoff        = (r2Unclamped < rCutoffMaxSq);

                    /* Pairs beyond the cut-off can only contribute via the perturbed-exclusion
                     * corrections, so included pairs beyond the cut-off can be skipped. */
                    if (pair.withinCutoff || !pair.pairIncluded)
                    {
                        // Ensure the distance does not become so small that r^-12 overflows
                        pair.r2    = sycl::max(r2Unclamped, c_nbnxmMinDistanceSquared);
                        pair.rInv  = sycl::rsqrt(pair.r2);
                        pair.r2Inv = pair.rInv * pair.rInv;

                        if (pair.pairIncluded && pair.withinCutoff)
                        {
                            nbfeComputeLjParams<useSoftCore, vdwType>(
                                    params, gm_nbfp, numTypes, typeiAB, typejAB, ljcpiAB, ljcpjAB, &pair);
                        }

                        NbfeAccumulators acc;
                        nbfeAccumulatePair<doCalcEnergies, useSoftCore, elecType, vdwType>(
                                params, lf, pair, gm_coulombTab, &acc);

                        if constexpr (doCalcEnergies)
                        {
                            totals.energyElec += acc.energyElec;
                            totals.energyVdw += acc.energyVdw;
                            totals.dvdlElec += acc.dvdlElec;
                            totals.dvdlVdw += acc.dvdlVdw;
                        }

                        if constexpr (!computeForeignLambda)
                        {
                            if (acc.scalarForcePerDistance != 0.0F)
                            {
                                const Float3 forceIJ = rv * acc.scalarForcePerDistance;
                                /* Accumulate the i-forces in registers */
                                fci += forceIJ;
                                /* Accumulate the j-forces directly in global memory */
                                staggeredAtomicAddForce(gm_f, -forceIJ, aj, laneIdx);
                            }
                        }
                    }
                } // if (j < nj1)
            } // for j-blocks

            if constexpr (!computeForeignLambda)
            {
                /* Reduce the i-forces over the sub-group and accumulate them in global memory */
                const Float3 fciSum = { sycl::reduce_over_group(sg, fci[XX], sycl::plus<float>()),
                                        sycl::reduce_over_group(sg, fci[YY], sycl::plus<float>()),
                                        sycl::reduce_over_group(sg, fci[ZZ], sycl::plus<float>()) };

                /* Skip central shifts when summing shift forces */
                const bool calcFShift = calcShift && (shiftIdx != gmx::c_centralShiftIndex);

                if (sg.leader())
                {
                    atomicFetchAdd(gm_f[ai], fciSum);
                    if (calcFShift)
                    {
                        atomicFetchAdd(gm_fShift[shiftIdx], fciSum);
                    }
                }
            }

            if constexpr (doCalcEnergies)
            {
                const float energyVdwSum =
                        sycl::reduce_over_group(sg, totals.energyVdw, sycl::plus<float>());
                const float energyElecSum =
                        sycl::reduce_over_group(sg, totals.energyElec, sycl::plus<float>());
                const float dvdlVdwSum =
                        sycl::reduce_over_group(sg, totals.dvdlVdw, sycl::plus<float>());
                const float dvdlElecSum =
                        sycl::reduce_over_group(sg, totals.dvdlElec, sycl::plus<float>());
                if (sg.leader())
                {
                    atomicFetchAdd(gm_energyVdw[lambdaIdx], energyVdwSum);
                    atomicFetchAdd(gm_energyElec[lambdaIdx], energyElecSum);
                    atomicFetchAdd(gm_dvdlVdw[lambdaIdx], dvdlVdwSum);
                    atomicFetchAdd(gm_dvdlElec[lambdaIdx], dvdlElecSum);
                }
            }
        }
    };
}

//! \brief Free-energy kernel launch code.
template<int subGroupSize, bool computeForeignLambda, bool doCalcEnergies, bool useSoftCore, enum ElecType elecType, enum VdwType vdwType, class... Args>
static void launchNbfeKernel(const DeviceStream& deviceStream, const int numiAtoms, Args&&... args)
{
    using kernelNameType =
            NbfeKernel<computeForeignLambda, doCalcEnergies, useSoftCore, elecType, vdwType, subGroupSize>;

    /* Kernel launch config: each sub-group handles one i-list entry. */
    constexpr int           workGroupSize = subGroupSize * sc_nbfeSubGroupsPerGroup;
    const int               numGroups     = divideRoundUp(numiAtoms, sc_nbfeSubGroupsPerGroup);
    const sycl::range<1>    blockSize{ workGroupSize };
    const sycl::range<1>    globalSize{ static_cast<size_t>(numGroups) * workGroupSize };
    const sycl::nd_range<1> range{ globalSize, blockSize };

    sycl::queue q = deviceStream.stream();

    auto kernelFunctionBuilder =
            nbfeKernel<subGroupSize, computeForeignLambda, doCalcEnergies, useSoftCore, elecType, vdwType>;
    syclSubmitWithoutCghOrEvent<kernelNameType>(
            q, kernelFunctionBuilder, range, std::forward<Args>(args)...);
}

//! \brief Select the templated free-energy kernel and launch it.
template<int subGroupSize, bool computeForeignLambda, bool doCalcEnergies, class... Args>
static void chooseAndLaunchNbfeKernel(enum ElecType elecType,
                                      enum VdwType  vdwType,
                                      const bool    useSoftCore,
                                      Args&&... args)
{
    gmx::dispatchTemplatedFunction(
            [&](auto elecType_, auto vdwType_, auto useSoftCore_)
            {
                return launchNbfeKernel<subGroupSize, computeForeignLambda, doCalcEnergies, useSoftCore_, elecType_, vdwType_>(
                        std::forward<Args>(args)...);
            },
            elecType,
            vdwType,
            useSoftCore);
}

template<int subGroupSize, bool computeForeignLambda, bool doCalcEnergies>
static void launchNbfeKernelHelperImpl(NbnxmGpu* nb, const bool calcShift, const InteractionLocality iloc)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    GpuFeplist*         feplist      = nb->feplist[iloc].get();
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    if (feplist->numiAtoms == 0)
    {
        /* Nothing to do for an empty list */
        return;
    }

    chooseAndLaunchNbfeKernel<subGroupSize, computeForeignLambda, doCalcEnergies>(
            nbp->elecType,
            nbp->vdwType,
            nbp->alphaVdw != 0.0F,
            deviceStream,
            feplist->numiAtoms,
            adat->xq.get_pointer(),
            adat->f.get_pointer(),
            adat->shiftVec.get_pointer(),
            adat->fShift.get_pointer(),
            (computeForeignLambda ? adat->eElecForeign : adat->eElec).get_pointer(),
            (computeForeignLambda ? adat->eLJForeign : adat->eLJ).get_pointer(),
            (computeForeignLambda ? adat->dvdlElecForeign : adat->dvdlElec).get_pointer(),
            (computeForeignLambda ? adat->dvdlLJForeign : adat->dvdlLJ).get_pointer(),
            adat->q4.get_pointer(),
            adat->atomTypes4.get_pointer(),
            adat->ljComb4.get_pointer(),
            nbp->nbfp.get_pointer(),
            nbp->coulomb_tab.get_pointer(),
            feplist->iinr.get_pointer(),
            feplist->jIndex.get_pointer(),
            feplist->jjnr.get_pointer(),
            feplist->shift.get_pointer(),
            feplist->exclFep.get_pointer(),
            feplist->numiAtoms,
            adat->numTypes,
            nbp->epsfac,
            nbp->lambdaPower,
            nbp->lambdaCoul,
            nbp->lambdaVdw,
            nbp->allLambdaCoul.get_pointer(),
            nbp->allLambdaVdw.get_pointer(),
            nb->fephostdata->allLambdaCoul.size(),
            makeNbfeInteractionParams(*nbp),
            calcShift);
}

template<int subGroupSize, bool doCalcEnergies>
void launchNbfeKernelHelper(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc)
{
    GMX_ASSERT(doCalcEnergies == stepWork.computeEnergy, "Wrong template called");
    launchNbfeKernelHelperImpl<subGroupSize, false, doCalcEnergies>(nb, stepWork.computeVirial, iloc);
}

template<int subGroupSize>
void launchNbfeForeignKernelHelper(NbnxmGpu* nb, InteractionLocality iloc)
{
    launchNbfeKernelHelperImpl<subGroupSize, true, true>(nb, false, iloc);
}

} // namespace gmx

#endif // GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_BODY_H
