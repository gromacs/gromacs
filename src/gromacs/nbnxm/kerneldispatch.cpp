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

#include "kerneldispatch.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
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

void
nbnxn_kernel_cpu(nonbonded_verlet_group_t  *nbvg,
                 nbnxn_atomdata_t          *nbat,
                 const interaction_const_t *ic,
                 rvec                      *shiftVectors,
                 int                        forceFlags,
                 int                        clearF,
                 real                      *fshift,
                 real                      *vCoulomb,
                 real                      *vVdw)
{

    int                      coulkt;
    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        coulkt = coulktRF;
    }
    else
    {
        if (nbvg->ewald_excl == ewaldexclTable)
        {
            if (ic->rcoulomb == ic->rvdw)
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
            if (ic->rcoulomb == ic->rvdw)
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
    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
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
    else if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == eljpmeGEOM)
        {
            vdwkt = vdwktLJEWALDCOMBGEOM;
        }
        else
        {
            vdwkt = vdwktLJEWALDCOMBLB;
            /* At setup we (should have) selected the C reference kernel */
            GMX_RELEASE_ASSERT(nbvg->kernel_type == nbnxnk4x4_PlainC, "Only the C reference nbnxn SIMD kernel supports LJ-PME with LB combination rules");
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unsupported VdW interaction type");
    }

    int                nnbl = nbvg->nbl_lists.nnbl;
    NbnxnPairlistCpu **nbl  = nbvg->nbl_lists.nbl;

    int gmx_unused     nthreads = gmx_omp_nthreads_get(emntNonbonded);
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
            switch (nbvg->kernel_type)
            {
                case nbnxnk4x4_PlainC:
                    nbnxn_kernel_noener_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                           ic,
                                                           shiftVectors,
                                                           out->f.data(),
                                                           fshift_p);
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case nbnxnk4xN_SIMD_2xNN:
                    nbnxm_kernel_noener_simd_2xmm[coulkt][vdwkt](nbl[nb], nbat,
                                                                 ic,
                                                                 shiftVectors,
                                                                 out->f.data(),
                                                                 fshift_p);
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case nbnxnk4xN_SIMD_4xN:
                    nbnxm_kernel_noener_simd_4xm[coulkt][vdwkt](nbl[nb], nbat,
                                                                ic,
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

            switch (nbvg->kernel_type)
            {
                case nbnxnk4x4_PlainC:
                    nbnxn_kernel_ener_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                         ic,
                                                         shiftVectors,
                                                         out->f.data(),
                                                         fshift_p,
                                                         out->Vvdw.data(),
                                                         out->Vc.data());
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case nbnxnk4xN_SIMD_2xNN:
                    nbnxm_kernel_ener_simd_2xmm[coulkt][vdwkt](nbl[nb], nbat,
                                                               ic,
                                                               shiftVectors,
                                                               out->f.data(),
                                                               fshift_p,
                                                               out->Vvdw.data(),
                                                               out->Vc.data());
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case nbnxnk4xN_SIMD_4xN:
                    nbnxm_kernel_ener_simd_4xm[coulkt][vdwkt](nbl[nb], nbat,
                                                              ic,
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

            switch (nbvg->kernel_type)
            {
                case nbnxnk4x4_PlainC:
                    unrollj = c_nbnxnCpuIClusterSize;
                    nbnxn_kernel_energrp_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                            ic,
                                                            shiftVectors,
                                                            out->f.data(),
                                                            fshift_p,
                                                            out->Vvdw.data(),
                                                            out->Vc.data());
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case nbnxnk4xN_SIMD_2xNN:
                    unrollj = GMX_SIMD_REAL_WIDTH/2;
                    nbnxm_kernel_energrp_simd_2xmm[coulkt][vdwkt](nbl[nb], nbat,
                                                                  ic,
                                                                  shiftVectors,
                                                                  out->f.data(),
                                                                  fshift_p,
                                                                  out->VSvdw.data(),
                                                                  out->VSc.data());
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case nbnxnk4xN_SIMD_4xN:
                    unrollj = GMX_SIMD_REAL_WIDTH;
                    nbnxm_kernel_energrp_simd_4xm[coulkt][vdwkt](nbl[nb], nbat,
                                                                 ic,
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

            if (nbvg->kernel_type != nbnxnk4x4_PlainC)
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
