/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#include "gmxpre.h"

#include "kernels_reference/kernel_gpu_ref.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "kernel_common.h"
#include "nbnxm_gpu.h"
#include "nbnxm_simd.h"
#include "pairlistset.h"
#include "pairlistsets.h"
#define INCLUDE_KERNELFUNCTION_TABLES
#include "kernels_reference/kernel_ref.h"
#if GMX_HAVE_NBNXM_SIMD_2XMM
#    include "kernels_simd_2xmm/kernels.h"
#endif
#if GMX_HAVE_NBNXM_SIMD_4XM
#    include "kernels_simd_4xm/kernels.h"
#endif
#undef INCLUDE_FUNCTION_TABLES

/*! \brief Clears the energy group output buffers
 *
 * \param[in,out] out  nbnxn kernel output struct
 */
static void clearGroupEnergies(nbnxn_atomdata_output_t* out)
{
    std::fill(out->Vvdw.begin(), out->Vvdw.end(), 0.0_real);
    std::fill(out->Vc.begin(), out->Vc.end(), 0.0_real);
    std::fill(out->VSvdw.begin(), out->VSvdw.end(), 0.0_real);
    std::fill(out->VSc.begin(), out->VSc.end(), 0.0_real);
}

/*! \brief Reduce the group-pair energy buffers produced by a SIMD kernel
 * to single terms in the output buffers.
 *
 * The SIMD kernels produce a large number of energy buffer in SIMD registers
 * to avoid scattered reads and writes.
 *
 * \tparam        unrollj         The unroll size for j-particles in the SIMD kernel
 * \param[in]     numGroups       The number of energy groups
 * \param[in]     numGroups_2log  Log2 of numGroups, rounded up
 * \param[in,out] out             Struct with energy buffers
 */
template<int unrollj>
static void reduceGroupEnergySimdBuffers(int numGroups, int numGroups_2log, nbnxn_atomdata_output_t* out)
{
    const int unrollj_half = unrollj / 2;
    /* Energies are stored in SIMD registers with size 2^numGroups_2log */
    const int numGroupsStorage = (1 << numGroups_2log);

    const real* gmx_restrict vVdwSimd     = out->VSvdw.data();
    const real* gmx_restrict vCoulombSimd = out->VSc.data();
    real* gmx_restrict       vVdw         = out->Vvdw.data();
    real* gmx_restrict       vCoulomb     = out->Vc.data();

    /* The size of the SIMD energy group buffer array is:
     * numGroups*numGroups*numGroupsStorage*unrollj_half*simd_width
     */
    for (int i = 0; i < numGroups; i++)
    {
        for (int j1 = 0; j1 < numGroups; j1++)
        {
            for (int j0 = 0; j0 < numGroups; j0++)
            {
                int c = ((i * numGroups + j1) * numGroupsStorage + j0) * unrollj_half * unrollj;
                for (int s = 0; s < unrollj_half; s++)
                {
                    vVdw[i * numGroups + j0] += vVdwSimd[c + 0];
                    vVdw[i * numGroups + j1] += vVdwSimd[c + 1];
                    vCoulomb[i * numGroups + j0] += vCoulombSimd[c + 0];
                    vCoulomb[i * numGroups + j1] += vCoulombSimd[c + 1];
                    c += unrollj + 2;
                }
            }
        }
    }
}

CoulombKernelType getCoulombKernelType(const Nbnxm::EwaldExclusionType ewaldExclusionType,
                                       const CoulombInteractionType    coulombInteractionType,
                                       const bool                      haveEqualCoulombVwdRadii)
{

    if (usingRF(coulombInteractionType) || coulombInteractionType == CoulombInteractionType::Cut)
    {
        return CoulombKernelType::ReactionField;
    }
    else
    {
        if (ewaldExclusionType == Nbnxm::EwaldExclusionType::Table)
        {
            if (haveEqualCoulombVwdRadii)
            {
                return CoulombKernelType::Table;
            }
            else
            {
                return CoulombKernelType::TableTwin;
            }
        }
        else
        {
            if (haveEqualCoulombVwdRadii)
            {
                return CoulombKernelType::Ewald;
            }
            else
            {
                return CoulombKernelType::EwaldTwin;
            }
        }
    }
}

int getVdwKernelType(const Nbnxm::KernelType    kernelType,
                     const LJCombinationRule    ljCombinationRule,
                     const VanDerWaalsType      vanDerWaalsType,
                     const InteractionModifiers interactionModifiers,
                     const LongRangeVdW         longRangeVdW)
{
    if (vanDerWaalsType == VanDerWaalsType::Cut)
    {
        switch (interactionModifiers)
        {
            case InteractionModifiers::None:
            case InteractionModifiers::PotShift:
                switch (ljCombinationRule)
                {
                    case LJCombinationRule::Geometric: return vdwktLJCUT_COMBGEOM;
                    case LJCombinationRule::LorentzBerthelot: return vdwktLJCUT_COMBLB;
                    case LJCombinationRule::None: return vdwktLJCUT_COMBNONE;
                    default: GMX_THROW(gmx::InvalidInputError("Unknown combination rule"));
                }
            case InteractionModifiers::ForceSwitch: return vdwktLJFORCESWITCH;
            case InteractionModifiers::PotSwitch: return vdwktLJPOTSWITCH;
            default:
                std::string errorMsg =
                        gmx::formatString("Unsupported VdW interaction modifier %s (%d)",
                                          enumValueToString(interactionModifiers),
                                          static_cast<int>(interactionModifiers));
                GMX_THROW(gmx::InvalidInputError(errorMsg));
        }
    }
    else if (vanDerWaalsType == VanDerWaalsType::Pme)
    {
        if (longRangeVdW == LongRangeVdW::Geom)
        {
            return vdwktLJEWALDCOMBGEOM;
        }
        else
        {
            /* At setup we (should have) selected the C reference kernel */
            GMX_RELEASE_ASSERT(kernelType == Nbnxm::KernelType::Cpu4x4_PlainC,
                               "Only the C reference nbnxn SIMD kernel supports LJ-PME with LB "
                               "combination rules");
            return vdwktLJEWALDCOMBLB;
        }
    }
    else
    {
        std::string errorMsg = gmx::formatString("Unsupported VdW interaction type %s (%d)",
                                                 enumValueToString(vanDerWaalsType),
                                                 static_cast<int>(vanDerWaalsType));
        GMX_THROW(gmx::InvalidInputError(errorMsg));
    }
}

/*! \brief Dispatches the non-bonded N versus M atom cluster CPU kernels.
 *
 * OpenMP parallelization is performed within this function.
 * Energy reduction, but not force and shift force reduction, is performed
 * within this function.
 *
 * \param[in]     pairlistSet   Pairlists with local or non-local interactions to compute
 * \param[in]     kernelSetup   The non-bonded kernel setup
 * \param[in,out] nbat          The atomdata for the interactions
 * \param[in]     ic            Non-bonded interaction constants
 * \param[in]     shiftVectors  The PBC shift vectors
 * \param[in]     stepWork      Flags that tell what to compute
 * \param[in]     clearF        Enum that tells if to clear the force output buffer
 * \param[out]    vCoulomb      Output buffer for Coulomb energies
 * \param[out]    vVdw          Output buffer for Van der Waals energies
 * \param[in]     wcycle        Pointer to cycle counting data structure.
 */
static void nbnxn_kernel_cpu(const PairlistSet&             pairlistSet,
                             const Nbnxm::KernelSetup&      kernelSetup,
                             nbnxn_atomdata_t*              nbat,
                             const interaction_const_t&     ic,
                             gmx::ArrayRef<const gmx::RVec> shiftVectors,
                             const gmx::StepWorkload&       stepWork,
                             int                            clearF,
                             real*                          vCoulomb,
                             real*                          vVdw,
                             gmx_wallcycle*                 wcycle)
{
    const nbnxn_atomdata_t::Params& nbatParams = nbat->params();

    GMX_ASSERT(ic.vdwtype != VanDerWaalsType::Pme
                       || ((ic.ljpme_comb_rule == LongRangeVdW::Geom
                            && nbatParams.ljCombinationRule == LJCombinationRule::Geometric)
                           || (ic.ljpme_comb_rule == LongRangeVdW::LB
                               && nbatParams.ljCombinationRule == LJCombinationRule::LorentzBerthelot)),
               "nbat combination rule parameters should match those for LJ-PME");

    const int coulkt = static_cast<int>(getCoulombKernelType(
            kernelSetup.ewaldExclusionType, ic.eeltype, (ic.rcoulomb == ic.rvdw)));
    const int vdwkt  = getVdwKernelType(
            kernelSetup.kernelType, nbatParams.ljCombinationRule, ic.vdwtype, ic.vdw_modifier, ic.ljpme_comb_rule);

    gmx::ArrayRef<const NbnxnPairlistCpu> pairlists = pairlistSet.cpuLists();

    const auto* shiftVecPointer = as_rvec_array(shiftVectors.data());

    int gmx_unused nthreads = gmx_omp_nthreads_get(ModuleMultiThread::Nonbonded);
    wallcycle_sub_start(wcycle, WallCycleSubCounter::NonbondedClear);
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (gmx::Index nb = 0; nb < pairlists.ssize(); nb++)
    {
        // Presently, the kernels do not call C++ code that can throw,
        // so no need for a try/catch pair in this OpenMP region.
        nbnxn_atomdata_output_t* out = &nbat->out[nb];

        if (clearF == enbvClearFYes)
        {
            clearForceBuffer(nbat, nb);

            clear_fshift(out->fshift.data());
        }

        if (nb == 0)
        {
            wallcycle_sub_stop(wcycle, WallCycleSubCounter::NonbondedClear);
            wallcycle_sub_start(wcycle, WallCycleSubCounter::NonbondedKernel);
        }

        // TODO: Change to reference
        const NbnxnPairlistCpu* pairlist = &pairlists[nb];

        if (!stepWork.computeEnergy)
        {
            /* Don't calculate energies */
            switch (kernelSetup.kernelType)
            {
                case Nbnxm::KernelType::Cpu4x4_PlainC:
                    nbnxn_kernel_noener_ref[coulkt][vdwkt](pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#if GMX_HAVE_NBNXM_SIMD_2XMM
                case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
                    gmx::nbnxmKernelNoenerSimd2xmm[coulkt][vdwkt](
                            pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#endif
#if GMX_HAVE_NBNXM_SIMD_4XM
                case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
                    gmx::nbnxmKernelNoenerSimd4xm[coulkt][vdwkt](pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#endif
                default: GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }
        }
        else if (out->Vvdw.size() == 1)
        {
            /* A single energy group (pair) */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            switch (kernelSetup.kernelType)
            {
                case Nbnxm::KernelType::Cpu4x4_PlainC:
                    nbnxn_kernel_ener_ref[coulkt][vdwkt](pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#if GMX_HAVE_NBNXM_SIMD_2XMM
                case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
                    gmx::nbnxmKernelEnerSimd2xmm[coulkt][vdwkt](pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#endif
#if GMX_HAVE_NBNXM_SIMD_4XM
                case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
                    gmx::nbnxmKernelEnerSimd4xm[coulkt][vdwkt](pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#endif
                default: GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }
        }
        else
        {
            /* Calculate energy group contributions */
            clearGroupEnergies(out);

            int unrollj = 0;

            switch (kernelSetup.kernelType)
            {
                case Nbnxm::KernelType::Cpu4x4_PlainC:
                    unrollj = c_nbnxnCpuIClusterSize;
                    nbnxn_kernel_energrp_ref[coulkt][vdwkt](pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#if GMX_HAVE_NBNXM_SIMD_2XMM
                case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
                    unrollj = GMX_SIMD_REAL_WIDTH / 2;
                    gmx::nbnxmKernelEnergrpSimd2xmm[coulkt][vdwkt](
                            pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#endif
#if GMX_HAVE_NBNXM_SIMD_4XM
                case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
                    unrollj = GMX_SIMD_REAL_WIDTH;
                    gmx::nbnxmKernelEnergrpSimd4xm[coulkt][vdwkt](
                            pairlist, nbat, &ic, shiftVecPointer, out);
                    break;
#endif
                default: GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }

            if (kernelSetup.kernelType != Nbnxm::KernelType::Cpu4x4_PlainC)
            {
                switch (unrollj)
                {
                    case 2:
                        reduceGroupEnergySimdBuffers<2>(nbatParams.nenergrp, nbatParams.neg_2log, out);
                        break;
                    case 4:
                        reduceGroupEnergySimdBuffers<4>(nbatParams.nenergrp, nbatParams.neg_2log, out);
                        break;
                    case 8:
                        reduceGroupEnergySimdBuffers<8>(nbatParams.nenergrp, nbatParams.neg_2log, out);
                        break;
                    default: GMX_RELEASE_ASSERT(false, "Unsupported j-unroll size");
                }
            }
        }
    }
    wallcycle_sub_stop(wcycle, WallCycleSubCounter::NonbondedKernel);

    if (stepWork.computeEnergy)
    {
        reduce_energies_over_lists(nbat, pairlists.ssize(), vVdw, vCoulomb);
    }
}

static void accountFlops(t_nrnb*                    nrnb,
                         const PairlistSet&         pairlistSet,
                         const nonbonded_verlet_t&  nbv,
                         const interaction_const_t& ic,
                         const gmx::StepWorkload&   stepWork)
{
    const bool usingGpuKernels = nbv.useGpu();

    int enr_nbnxn_kernel_ljc = eNRNB;
    if (usingRF(ic.eeltype) || ic.eeltype == CoulombInteractionType::Cut)
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_RF;
    }
    else if ((!usingGpuKernels && nbv.kernelSetup().ewaldExclusionType == Nbnxm::EwaldExclusionType::Analytical)
             || (usingGpuKernels && Nbnxm::gpu_is_kernel_ewald_analytical(nbv.gpuNbv())))
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_EWALD;
    }
    else
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_TAB;
    }
    int enr_nbnxn_kernel_lj = eNR_NBNXN_LJ;
    if (stepWork.computeEnergy)
    {
        /* In eNR_??? the nbnxn F+E kernels are always the F kernel + 1 */
        enr_nbnxn_kernel_ljc += 1;
        enr_nbnxn_kernel_lj += 1;
    }

    inc_nrnb(nrnb, enr_nbnxn_kernel_ljc, pairlistSet.natpair_ljq_);
    inc_nrnb(nrnb, enr_nbnxn_kernel_lj, pairlistSet.natpair_lj_);
    /* The Coulomb-only kernels are offset -eNR_NBNXN_LJ_RF+eNR_NBNXN_RF */
    inc_nrnb(nrnb, enr_nbnxn_kernel_ljc - eNR_NBNXN_LJ_RF + eNR_NBNXN_RF, pairlistSet.natpair_q_);

    if (ic.vdw_modifier == InteractionModifiers::ForceSwitch)
    {
        /* We add up the switch cost separately */
        inc_nrnb(nrnb,
                 eNR_NBNXN_ADD_LJ_FSW + (stepWork.computeEnergy ? 1 : 0),
                 pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_);
    }
    if (ic.vdw_modifier == InteractionModifiers::PotSwitch)
    {
        /* We add up the switch cost separately */
        inc_nrnb(nrnb,
                 eNR_NBNXN_ADD_LJ_PSW + (stepWork.computeEnergy ? 1 : 0),
                 pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_);
    }
    if (ic.vdwtype == VanDerWaalsType::Pme)
    {
        /* We add up the LJ Ewald cost separately */
        inc_nrnb(nrnb,
                 eNR_NBNXN_ADD_LJ_EWALD + (stepWork.computeEnergy ? 1 : 0),
                 pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_);
    }
}

void nonbonded_verlet_t::dispatchNonbondedKernel(gmx::InteractionLocality       iLocality,
                                                 const interaction_const_t&     ic,
                                                 const gmx::StepWorkload&       stepWork,
                                                 int                            clearF,
                                                 gmx::ArrayRef<const gmx::RVec> shiftvec,
                                                 gmx::ArrayRef<real> repulsionDispersionSR,
                                                 gmx::ArrayRef<real> CoulombSR,
                                                 t_nrnb*             nrnb) const
{
    const PairlistSet& pairlistSet = pairlistSets().pairlistSet(iLocality);

    switch (kernelSetup().kernelType)
    {
        case Nbnxm::KernelType::Cpu4x4_PlainC:
        case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
        case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
            nbnxn_kernel_cpu(pairlistSet,
                             kernelSetup(),
                             nbat_.get(),
                             ic,
                             shiftvec,
                             stepWork,
                             clearF,
                             CoulombSR.data(),
                             repulsionDispersionSR.data(),
                             wcycle_);
            break;

        case Nbnxm::KernelType::Gpu8x8x8:
            Nbnxm::gpu_launch_kernel(gpuNbv_, stepWork, iLocality);
            break;

        case Nbnxm::KernelType::Cpu8x8x8_PlainC:
            nbnxn_kernel_gpu_ref(pairlistSet.gpuList(),
                                 nbat_.get(),
                                 &ic,
                                 shiftvec,
                                 stepWork,
                                 clearF,
                                 nbat_->out[0].f,
                                 nbat_->out[0].fshift.data(),
                                 CoulombSR.data(),
                                 repulsionDispersionSR.data());
            break;

        default: GMX_RELEASE_ASSERT(false, "Invalid nonbonded kernel type passed!");
    }

    if (nrnb)
    {
        accountFlops(nrnb, pairlistSet, *this, ic, stepWork);
    }
}
