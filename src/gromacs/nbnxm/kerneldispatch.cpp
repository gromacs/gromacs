/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#include "gmxpre.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/kernels_reference/kernel_gpu_ref.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "kernel_common.h"
#define INCLUDE_KERNELFUNCTION_TABLES
#include "gromacs/nbnxm/kernels_reference/kernel_ref.h"
#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/nbnxm/kernels_simd_2xmm/kernels.h"
#endif
#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/nbnxm/kernels_simd_4xm/kernels.h"
#endif
#undef INCLUDE_FUNCTION_TABLES

/*! \brief Clears the energy group output buffers
 *
 * \param[in,out] out  nbnxn kernel output struct
 */
static void clearGroupEnergies(nbnxn_atomdata_output_t *out)
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
template <int unrollj> static void
reduceGroupEnergySimdBuffers(int                       numGroups,
                             int                       numGroups_2log,
                             nbnxn_atomdata_output_t  *out)
{
    const int                 unrollj_half     = unrollj/2;
    /* Energies are stored in SIMD registers with size 2^numGroups_2log */
    const int                 numGroupsStorage = (1 << numGroups_2log);

    const real * gmx_restrict vVdwSimd     = out->VSvdw.data();
    const real * gmx_restrict vCoulombSimd = out->VSc.data();
    real * gmx_restrict       vVdw         = out->Vvdw.data();
    real * gmx_restrict       vCoulomb     = out->Vc.data();

    /* The size of the SIMD energy group buffer array is:
     * numGroups*numGroups*numGroupsStorage*unrollj_half*simd_width
     */
    for (int i = 0; i < numGroups; i++)
    {
        for (int j1 = 0; j1 < numGroups; j1++)
        {
            for (int j0 = 0; j0 < numGroups; j0++)
            {
                int c = ((i*numGroups + j1)*numGroupsStorage + j0)*unrollj_half*unrollj;
                for (int s = 0; s < unrollj_half; s++)
                {
                    vVdw    [i*numGroups + j0] += vVdwSimd    [c + 0];
                    vVdw    [i*numGroups + j1] += vVdwSimd    [c + 1];
                    vCoulomb[i*numGroups + j0] += vCoulombSimd[c + 0];
                    vCoulomb[i*numGroups + j1] += vCoulombSimd[c + 1];
                    c                          += unrollj + 2;
                }
            }
        }
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
 * \param[in]     forceFlags    Flags that tell what to compute
 * \param[in]     clearF        Enum that tells if to clear the force output buffer
 * \param[out]    fshift        Shift force output buffer
 * \param[out]    vCoulomb      Output buffer for Coulomb energies
 * \param[out]    vVdw          Output buffer for Van der Waals energies
 */
static void
nbnxn_kernel_cpu(const nbnxn_pairlist_set_t     &pairlistSet,
                 const Nbnxm::KernelSetup       &kernelSetup,
                 nbnxn_atomdata_t               *nbat,
                 const interaction_const_t      &ic,
                 rvec                           *shiftVectors,
                 int                             forceFlags,
                 int                             clearF,
                 real                           *fshift,
                 real                           *vCoulomb,
                 real                           *vVdw)
{

    int                      coulkt;
    if (EEL_RF(ic.eeltype) || ic.eeltype == eelCUT)
    {
        coulkt = coulktRF;
    }
    else
    {
        if (kernelSetup.ewaldExclusionType == Nbnxm::EwaldExclusionType::Table)
        {
            if (ic.rcoulomb == ic.rvdw)
            {
                coulkt = coulktTAB;
            }
            else
            {
                coulkt = coulktTAB_TWIN;
            }
        }
        else
        {
            if (ic.rcoulomb == ic.rvdw)
            {
                coulkt = coulktEWALD;
            }
            else
            {
                coulkt = coulktEWALD_TWIN;
            }
        }
    }

    const nbnxn_atomdata_t::Params &nbatParams = nbat->params();

    int vdwkt = 0;
    if (ic.vdwtype == evdwCUT)
    {
        switch (ic.vdw_modifier)
        {
            case eintmodNONE:
            case eintmodPOTSHIFT:
                switch (nbatParams.comb_rule)
                {
                    case ljcrGEOM: vdwkt = vdwktLJCUT_COMBGEOM; break;
                    case ljcrLB:   vdwkt = vdwktLJCUT_COMBLB;   break;
                    case ljcrNONE: vdwkt = vdwktLJCUT_COMBNONE; break;
                    default:
                        GMX_RELEASE_ASSERT(false, "Unknown combination rule");
                }
                break;
            case eintmodFORCESWITCH:
                vdwkt = vdwktLJFORCESWITCH;
                break;
            case eintmodPOTSWITCH:
                vdwkt = vdwktLJPOTSWITCH;
                break;
            default:
                GMX_RELEASE_ASSERT(false, "Unsupported VdW interaction modifier");
        }
    }
    else if (ic.vdwtype == evdwPME)
    {
        if (ic.ljpme_comb_rule == eljpmeGEOM)
        {
            vdwkt = vdwktLJEWALDCOMBGEOM;
        }
        else
        {
            vdwkt = vdwktLJEWALDCOMBLB;
            /* At setup we (should have) selected the C reference kernel */
            GMX_RELEASE_ASSERT(kernelSetup.kernelType == Nbnxm::KernelType::Cpu4x4_PlainC, "Only the C reference nbnxn SIMD kernel supports LJ-PME with LB combination rules");
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unsupported VdW interaction type");
    }

    int                        nnbl = pairlistSet.nnbl;
    NbnxnPairlistCpu * const * nbl  = pairlistSet.nbl;

    int gmx_unused             nthreads = gmx_omp_nthreads_get(emntNonbonded);
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int nb = 0; nb < nnbl; nb++)
    {
        // Presently, the kernels do not call C++ code that can throw,
        // so no need for a try/catch pair in this OpenMP region.
        nbnxn_atomdata_output_t *out = &nbat->out[nb];

        if (clearF == enbvClearFYes)
        {
            clear_f(nbat, nb, out->f.data());
        }

        real *fshift_p;
        if ((forceFlags & GMX_FORCE_VIRIAL) && nnbl == 1)
        {
            fshift_p = fshift;
        }
        else
        {
            fshift_p = out->fshift.data();

            if (clearF == enbvClearFYes)
            {
                clear_fshift(fshift_p);
            }
        }

        if (!(forceFlags & GMX_FORCE_ENERGY))
        {
            /* Don't calculate energies */
            switch (kernelSetup.kernelType)
            {
                case Nbnxm::KernelType::Cpu4x4_PlainC:
                    nbnxn_kernel_noener_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                           &ic,
                                                           shiftVectors,
                                                           out->f.data(),
                                                           fshift_p);
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
                    nbnxm_kernel_noener_simd_2xmm[coulkt][vdwkt](nbl[nb], nbat,
                                                                 &ic,
                                                                 shiftVectors,
                                                                 out->f.data(),
                                                                 fshift_p);
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
                    nbnxm_kernel_noener_simd_4xm[coulkt][vdwkt](nbl[nb], nbat,
                                                                &ic,
                                                                shiftVectors,
                                                                out->f.data(),
                                                                fshift_p);
                    break;
#endif
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
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
                    nbnxn_kernel_ener_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                         &ic,
                                                         shiftVectors,
                                                         out->f.data(),
                                                         fshift_p,
                                                         out->Vvdw.data(),
                                                         out->Vc.data());
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
                    nbnxm_kernel_ener_simd_2xmm[coulkt][vdwkt](nbl[nb], nbat,
                                                               &ic,
                                                               shiftVectors,
                                                               out->f.data(),
                                                               fshift_p,
                                                               out->Vvdw.data(),
                                                               out->Vc.data());
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
                    nbnxm_kernel_ener_simd_4xm[coulkt][vdwkt](nbl[nb], nbat,
                                                              &ic,
                                                              shiftVectors,
                                                              out->f.data(),
                                                              fshift_p,
                                                              out->Vvdw.data(),
                                                              out->Vc.data());
                    break;
#endif
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
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
                    nbnxn_kernel_energrp_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                            &ic,
                                                            shiftVectors,
                                                            out->f.data(),
                                                            fshift_p,
                                                            out->Vvdw.data(),
                                                            out->Vc.data());
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
                    unrollj = GMX_SIMD_REAL_WIDTH/2;
                    nbnxm_kernel_energrp_simd_2xmm[coulkt][vdwkt](nbl[nb], nbat,
                                                                  &ic,
                                                                  shiftVectors,
                                                                  out->f.data(),
                                                                  fshift_p,
                                                                  out->VSvdw.data(),
                                                                  out->VSc.data());
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
                    unrollj = GMX_SIMD_REAL_WIDTH;
                    nbnxm_kernel_energrp_simd_4xm[coulkt][vdwkt](nbl[nb], nbat,
                                                                 &ic,
                                                                 shiftVectors,
                                                                 out->f.data(),
                                                                 fshift_p,
                                                                 out->VSvdw.data(),
                                                                 out->VSc.data());
                    break;
#endif
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }

            if (kernelSetup.kernelType != Nbnxm::KernelType::Cpu4x4_PlainC)
            {
                switch (unrollj)
                {
                    case 2:
                        reduceGroupEnergySimdBuffers<2>(nbatParams.nenergrp,
                                                        nbatParams.neg_2log,
                                                        out);
                        break;
                    case 4:
                        reduceGroupEnergySimdBuffers<4>(nbatParams.nenergrp,
                                                        nbatParams.neg_2log,
                                                        out);
                        break;
                    case 8:
                        reduceGroupEnergySimdBuffers<8>(nbatParams.nenergrp,
                                                        nbatParams.neg_2log,
                                                        out);
                        break;
                    default:
                        GMX_RELEASE_ASSERT(false, "Unsupported j-unroll size");
                }
            }
        }
    }

    if (forceFlags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(nbat, nnbl, vVdw, vCoulomb);
    }
}

static void accountFlops(t_nrnb                           *nrnb,
                         const nbnxn_pairlist_set_t       &pairlistSet,
                         const nonbonded_verlet_t         &nbv,
                         const interaction_const_t        &ic,
                         const int                         forceFlags)
{
    const bool usingGpuKernels = nbv.useGpu();

    int        enr_nbnxn_kernel_ljc;
    if (EEL_RF(ic.eeltype) || ic.eeltype == eelCUT)
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_RF;
    }
    else if ((!usingGpuKernels && nbv.kernelSetup().ewaldExclusionType == Nbnxm::EwaldExclusionType::Analytical) ||
             (usingGpuKernels && Nbnxm::gpu_is_kernel_ewald_analytical(nbv.gpu_nbv)))
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_EWALD;
    }
    else
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_TAB;
    }
    int enr_nbnxn_kernel_lj = eNR_NBNXN_LJ;
    if (forceFlags & GMX_FORCE_ENERGY)
    {
        /* In eNR_??? the nbnxn F+E kernels are always the F kernel + 1 */
        enr_nbnxn_kernel_ljc += 1;
        enr_nbnxn_kernel_lj  += 1;
    }

    inc_nrnb(nrnb, enr_nbnxn_kernel_ljc,
             pairlistSet.natpair_ljq);
    inc_nrnb(nrnb, enr_nbnxn_kernel_lj,
             pairlistSet.natpair_lj);
    /* The Coulomb-only kernels are offset -eNR_NBNXN_LJ_RF+eNR_NBNXN_RF */
    inc_nrnb(nrnb, enr_nbnxn_kernel_ljc-eNR_NBNXN_LJ_RF+eNR_NBNXN_RF,
             pairlistSet.natpair_q);

    const bool calcEnergy = ((forceFlags & GMX_FORCE_ENERGY) != 0);
    if (ic.vdw_modifier == eintmodFORCESWITCH)
    {
        /* We add up the switch cost separately */
        inc_nrnb(nrnb, eNR_NBNXN_ADD_LJ_FSW + (calcEnergy ? 1 : 0),
                 pairlistSet.natpair_ljq + pairlistSet.natpair_lj);
    }
    if (ic.vdw_modifier == eintmodPOTSWITCH)
    {
        /* We add up the switch cost separately */
        inc_nrnb(nrnb, eNR_NBNXN_ADD_LJ_PSW + (calcEnergy ? 1 : 0),
                 pairlistSet.natpair_ljq + pairlistSet.natpair_lj);
    }
    if (ic.vdwtype == evdwPME)
    {
        /* We add up the LJ Ewald cost separately */
        inc_nrnb(nrnb, eNR_NBNXN_ADD_LJ_EWALD + (calcEnergy ? 1 : 0),
                 pairlistSet.natpair_ljq + pairlistSet.natpair_lj);
    }
}

void
nonbonded_verlet_t::dispatchNonbondedKernel(Nbnxm::InteractionLocality iLocality,
                                            const interaction_const_t &ic,
                                            int                        forceFlags,
                                            int                        clearF,
                                            t_forcerec                *fr,
                                            gmx_enerdata_t            *enerd,
                                            t_nrnb                    *nrnb)
{
    const nbnxn_pairlist_set_t &pairlistSet = pairlistSets().pairlistSet(iLocality);

    switch (kernelSetup().kernelType)
    {
        case Nbnxm::KernelType::Cpu4x4_PlainC:
        case Nbnxm::KernelType::Cpu4xN_Simd_4xN:
        case Nbnxm::KernelType::Cpu4xN_Simd_2xNN:
            nbnxn_kernel_cpu(pairlistSet,
                             kernelSetup(),
                             nbat.get(),
                             ic,
                             fr->shift_vec,
                             forceFlags,
                             clearF,
                             fr->fshift[0],
                             enerd->grpp.ener[egCOULSR],
                             fr->bBHAM ?
                             enerd->grpp.ener[egBHAMSR] :
                             enerd->grpp.ener[egLJSR]);
            break;

        case Nbnxm::KernelType::Gpu8x8x8:
            Nbnxm::gpu_launch_kernel(gpu_nbv, forceFlags, iLocality);
            break;

        case Nbnxm::KernelType::Cpu8x8x8_PlainC:
            nbnxn_kernel_gpu_ref(pairlistSet.nblGpu[0],
                                 nbat.get(), &ic,
                                 fr->shift_vec,
                                 forceFlags,
                                 clearF,
                                 nbat->out[0].f,
                                 fr->fshift[0],
                                 enerd->grpp.ener[egCOULSR],
                                 fr->bBHAM ?
                                 enerd->grpp.ener[egBHAMSR] :
                                 enerd->grpp.ener[egLJSR]);
            break;

        default:
            GMX_RELEASE_ASSERT(false, "Invalid nonbonded kernel type passed!");

    }

    accountFlops(nrnb, pairlistSet, *this, ic, forceFlags);
}

void
nonbonded_verlet_t::dispatchFreeEnergyKernel(Nbnxm::InteractionLocality  iLocality,
                                             t_forcerec                 *fr,
                                             rvec                        x[],
                                             rvec                        f[],
                                             const t_mdatoms            &mdatoms,
                                             t_lambda                   *fepvals,
                                             real                       *lambda,
                                             gmx_enerdata_t             *enerd,
                                             const int                   forceFlags,
                                             t_nrnb                     *nrnb)
{
    const gmx::ArrayRef<t_nblist const * const > nbl_fep = pairlistSets().pairlistSet(iLocality).nbl_fep;

    /* When the first list is empty, all are empty and there is nothing to do */
    if (nbl_fep[0]->nrj == 0)
    {
        return;
    }

    int donb_flags = 0;
    /* Add short-range interactions */
    donb_flags |= GMX_NONBONDED_DO_SR;

    /* Currently all group scheme kernels always calculate (shift-)forces */
    if (forceFlags & GMX_FORCE_FORCES)
    {
        donb_flags |= GMX_NONBONDED_DO_FORCE;
    }
    if (forceFlags & GMX_FORCE_VIRIAL)
    {
        donb_flags |= GMX_NONBONDED_DO_SHIFTFORCE;
    }
    if (forceFlags & GMX_FORCE_ENERGY)
    {
        donb_flags |= GMX_NONBONDED_DO_POTENTIAL;
    }

    nb_kernel_data_t kernel_data;
    real             dvdl_nb[efptNR] = { 0 };
    kernel_data.flags  = donb_flags;
    kernel_data.lambda = lambda;
    kernel_data.dvdl   = dvdl_nb;

    kernel_data.energygrp_elec = enerd->grpp.ener[egCOULSR];
    kernel_data.energygrp_vdw  = enerd->grpp.ener[egLJSR];

    GMX_ASSERT(gmx_omp_nthreads_get(emntNonbonded) == nbl_fep.ssize(), "Number of lists should be same as number of NB threads");

#pragma omp parallel for schedule(static) num_threads(nbl_fep.ssize())
    for (int th = 0; th < nbl_fep.ssize(); th++)
    {
        try
        {
            gmx_nb_free_energy_kernel(nbl_fep[th],
                                      x, f, fr, &mdatoms, &kernel_data, nrnb);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    if (fepvals->sc_alpha != 0)
    {
        enerd->dvdl_nonlin[efptVDW]  += dvdl_nb[efptVDW];
        enerd->dvdl_nonlin[efptCOUL] += dvdl_nb[efptCOUL];
    }
    else
    {
        enerd->dvdl_lin[efptVDW]  += dvdl_nb[efptVDW];
        enerd->dvdl_lin[efptCOUL] += dvdl_nb[efptCOUL];
    }

    /* If we do foreign lambda and we have soft-core interactions
     * we have to recalculate the (non-linear) energies contributions.
     */
    if (fepvals->n_lambda > 0 && (forceFlags & GMX_FORCE_DHDL) && fepvals->sc_alpha != 0)
    {
        real lam_i[efptNR];
        kernel_data.flags          = (donb_flags & ~(GMX_NONBONDED_DO_FORCE | GMX_NONBONDED_DO_SHIFTFORCE)) | GMX_NONBONDED_DO_FOREIGNLAMBDA;
        kernel_data.lambda         = lam_i;
        kernel_data.energygrp_elec = enerd->foreign_grpp.ener[egCOULSR];
        kernel_data.energygrp_vdw  = enerd->foreign_grpp.ener[egLJSR];
        /* Note that we add to kernel_data.dvdl, but ignore the result */

        for (int i = 0; i < enerd->n_lambda; i++)
        {
            for (int j = 0; j < efptNR; j++)
            {
                lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
            }
            reset_foreign_enerdata(enerd);
#pragma omp parallel for schedule(static) num_threads(nbl_fep.ssize())
            for (int th = 0; th < nbl_fep.ssize(); th++)
            {
                try
                {
                    gmx_nb_free_energy_kernel(nbl_fep[th],
                                              x, f, fr, &mdatoms, &kernel_data, nrnb);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

            sum_epot(&(enerd->foreign_grpp), enerd->foreign_term);
            enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
        }
    }
}
