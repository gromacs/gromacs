/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "config.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
#define INCLUDE_KERNELFUNCTION_TABLES
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref.h"
#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn.h"
#endif
#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn.h"
#endif
#undef INCLUDE_FUNCTION_TABLES
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"


template <int unrollj> static void
reduce_group_energies(int ng, int ng_2log,
                      const real *VSvdw, const real *VSc,
                      real * gmx_restrict Vvdw,
                      real * gmx_restrict Vc)
{
    const int unrollj_half = unrollj/2;

    const int ng_p2        = (1 << ng_2log);

    /* The size of the x86 SIMD energy group buffer array is:
     * ng*ng*ng_p2*unrollj_half*simd_width
     */
    for (int i = 0; i < ng; i++)
    {
        for (int j = 0; j < ng; j++)
        {
            Vvdw[i*ng + j] = 0;
            Vc[i*ng + j]   = 0;
        }

        for (int j1 = 0; j1 < ng; j1++)
        {
            for (int j0 = 0; j0 < ng; j0++)
            {
                int c = ((i*ng + j1)*ng_p2 + j0)*unrollj_half*unrollj;
                for (int s = 0; s < unrollj_half; s++)
                {
                    Vvdw[i*ng+j0] += VSvdw[c+0];
                    Vvdw[i*ng+j1] += VSvdw[c+1];
                    Vc  [i*ng+j0] += VSc  [c+0];
                    Vc  [i*ng+j1] += VSc  [c+1];
                    c             += unrollj + 2;
                }
            }
        }
    }
}

void
nbnxn_kernel(nonbonded_verlet_group_t  *nbvg,
             const interaction_const_t *ic,
             rvec                      *shift_vec,
             int                        force_flags,
             int                        clearF,
             real                      *fshift,
             real                      *Vc,
             real                      *Vvdw)
{
    int                      nnbl = nbvg->nbl_lists.nnbl;
    nbnxn_pairlist_t       **nbl  = nbvg->nbl_lists.nbl;
    const nbnxn_atomdata_t  *nbat = nbvg->nbat;

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

    int vdwkt = 0;
    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodNONE:
            case eintmodPOTSHIFT:
                switch (nbat->comb_rule)
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
            /* At setup we selected the C reference kernel */
            GMX_RELEASE_ASSERT(nbvg->kernel_type == nbnxnk4x4_PlainC, "Only the C reference nbnxn SIMD kernel supports LJ-PME with LB combination rules");
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unsupported VdW interaction type");
    }
    // cppcheck-suppress unreadVariable
    int nthreads = gmx_omp_nthreads_get(emntNonbonded);
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int nb = 0; nb < nnbl; nb++)
    {
        // Presently, the kernels do not call C++ code that can throw, so
        // no need for a try/catch pair in this OpenMP region.
        nbnxn_atomdata_output_t *out;
        real                    *fshift_p;

        out = &nbat->out[nb];

        if (clearF == enbvClearFYes)
        {
            clear_f(nbat, nb, out->f);
        }

        if ((force_flags & GMX_FORCE_VIRIAL) && nnbl == 1)
        {
            fshift_p = fshift;
        }
        else
        {
            fshift_p = out->fshift;

            if (clearF == enbvClearFYes)
            {
                clear_fshift(fshift_p);
            }
        }

        if (!(force_flags & GMX_FORCE_ENERGY))
        {
            /* Don't calculate energies */
            switch (nbvg->kernel_type)
            {
                case nbnxnk4x4_PlainC:
                    nbnxn_kernel_noener_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                           ic,
                                                           shift_vec,
                                                           out->f,
                                                           fshift_p);
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case nbnxnk4xN_SIMD_2xNN:
                    nbnxn_kernel_noener_simd_2xnn[coulkt][vdwkt](nbl[nb], nbat,
                                                                 ic,
                                                                 shift_vec,
                                                                 out->f,
                                                                 fshift_p);
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case nbnxnk4xN_SIMD_4xN:
                    nbnxn_kernel_noener_simd_4xn[coulkt][vdwkt](nbl[nb], nbat,
                                                                ic,
                                                                shift_vec,
                                                                out->f,
                                                                fshift_p);
                    break;
#endif
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }
        }
        else if (out->nV == 1)
        {
            /* No energy groups */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            switch (nbvg->kernel_type)
            {
                case nbnxnk4x4_PlainC:
                    nbnxn_kernel_ener_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                         ic,
                                                         shift_vec,
                                                         out->f,
                                                         fshift_p,
                                                         out->Vvdw,
                                                         out->Vc);
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case nbnxnk4xN_SIMD_2xNN:
                    nbnxn_kernel_ener_simd_2xnn[coulkt][vdwkt](nbl[nb], nbat,
                                                               ic,
                                                               shift_vec,
                                                               out->f,
                                                               fshift_p,
                                                               out->Vvdw,
                                                               out->Vc);
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case nbnxnk4xN_SIMD_4xN:
                    nbnxn_kernel_ener_simd_4xn[coulkt][vdwkt](nbl[nb], nbat,
                                                              ic,
                                                              shift_vec,
                                                              out->f,
                                                              fshift_p,
                                                              out->Vvdw,
                                                              out->Vc);
                    break;
#endif
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }
        }
        else
        {
            /* Calculate energy group contributions */
            int i;

            for (i = 0; i < out->nVS; i++)
            {
                out->VSvdw[i] = 0;
            }
            for (i = 0; i < out->nVS; i++)
            {
                out->VSc[i] = 0;
            }

            int unrollj = 0;

            switch (nbvg->kernel_type)
            {
                case nbnxnk4x4_PlainC:
                    unrollj = NBNXN_CPU_CLUSTER_I_SIZE;
                    nbnxn_kernel_energrp_ref[coulkt][vdwkt](nbl[nb], nbat,
                                                            ic,
                                                            shift_vec,
                                                            out->f,
                                                            fshift_p,
                                                            out->VSvdw,
                                                            out->VSc);
                    break;
#ifdef GMX_NBNXN_SIMD_2XNN
                case nbnxnk4xN_SIMD_2xNN:
                    unrollj = GMX_SIMD_REAL_WIDTH/2;
                    nbnxn_kernel_energrp_simd_2xnn[coulkt][vdwkt](nbl[nb], nbat,
                                                                  ic,
                                                                  shift_vec,
                                                                  out->f,
                                                                  fshift_p,
                                                                  out->VSvdw,
                                                                  out->VSc);
                    break;
#endif
#ifdef GMX_NBNXN_SIMD_4XN
                case nbnxnk4xN_SIMD_4xN:
                    unrollj = GMX_SIMD_REAL_WIDTH;
                    nbnxn_kernel_energrp_simd_4xn[coulkt][vdwkt](nbl[nb], nbat,
                                                                 ic,
                                                                 shift_vec,
                                                                 out->f,
                                                                 fshift_p,
                                                                 out->VSvdw,
                                                                 out->VSc);
                    break;
#endif
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported kernel architecture");
            }

            switch (unrollj)
            {
                case 2:
                    reduce_group_energies<2>(nbat->nenergrp, nbat->neg_2log,
                                             out->VSvdw, out->VSc,
                                             out->Vvdw, out->Vc);
                    break;
                case 4:
                    reduce_group_energies<4>(nbat->nenergrp, nbat->neg_2log,
                                             out->VSvdw, out->VSc,
                                             out->Vvdw, out->Vc);
                    break;
                case 8:
                    reduce_group_energies<8>(nbat->nenergrp, nbat->neg_2log,
                                             out->VSvdw, out->VSc,
                                             out->Vvdw, out->Vc);
                    break;
                default:
                    GMX_RELEASE_ASSERT(false, "Unsupported j-unroll size");
            }
        }
    }

    if (force_flags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(nbat, nnbl, Vvdw, Vc);
    }
}
