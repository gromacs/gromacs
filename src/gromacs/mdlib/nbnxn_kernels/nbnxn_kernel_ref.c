/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "typedefs.h"
#include "vec.h"
#include "smalloc.h"
#include "force.h"
#include "gmx_omp_nthreads.h"
#include "nbnxn_kernel_ref.h"
#include "../nbnxn_consts.h"
#include "nbnxn_kernel_common.h"

/*! \brief Typedefs for declaring lookup tables of kernel functions.
 */

typedef void (*p_nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                                  const nbnxn_atomdata_t     *nbat,
                                  const interaction_const_t  *ic,
                                  rvec                       *shift_vec,
                                  real                       *f,
                                  real                       *fshift);

typedef void (*p_nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                                const nbnxn_atomdata_t     *nbat,
                                const interaction_const_t  *ic,
                                rvec                       *shift_vec,
                                real                       *f,
                                real                       *fshift,
                                real                       *Vvdw,
                                real                       *Vc);

/* Analytical reaction-field kernels */
#define CALC_COUL_RF
#include "nbnxn_kernel_ref_includes.h"
#define VDW_FORCE_SWITCH
#include "nbnxn_kernel_ref_includes.h"
#undef VDW_FORCE_SWITCH
#define VDW_POT_SWITCH
#include "nbnxn_kernel_ref_includes.h"
#undef VDW_POT_SWITCH
#undef CALC_COUL_RF


/* Tabulated exclusion interaction electrostatics kernels */
#define CALC_COUL_TAB
#include "nbnxn_kernel_ref_includes.h"
#define VDW_FORCE_SWITCH
#include "nbnxn_kernel_ref_includes.h"
#undef VDW_FORCE_SWITCH
#define VDW_POT_SWITCH
#include "nbnxn_kernel_ref_includes.h"
#undef VDW_POT_SWITCH
/* Twin-range cut-off kernels */
#define VDW_CUTOFF_CHECK
#include "nbnxn_kernel_ref_includes.h"
#define VDW_FORCE_SWITCH
#include "nbnxn_kernel_ref_includes.h"
#undef VDW_FORCE_SWITCH
#define VDW_POT_SWITCH
#include "nbnxn_kernel_ref_includes.h"
#undef VDW_POT_SWITCH
#undef VDW_CUTOFF_CHECK
#undef CALC_COUL_TAB


enum {
    coultRF, coultTAB, coultTAB_TWIN, coultNR
};

enum {
    vdwtCUT, vdwtFSWITCH, vdwtPSWITCH, vdwtNR
};

p_nbk_func_noener p_nbk_c_noener[coultNR][vdwtNR] =
{
    { nbnxn_kernel_ref_rf_lj_cut_noener,       nbnxn_kernel_ref_rf_lj_fswitch_noener,       nbnxn_kernel_ref_rf_lj_pswitch_noener },
    { nbnxn_kernel_ref_tab_lj_cut_noener,      nbnxn_kernel_ref_tab_lj_fswitch_noener,      nbnxn_kernel_ref_tab_lj_pswitch_noener },
    { nbnxn_kernel_ref_tab_twin_lj_cut_noener, nbnxn_kernel_ref_tab_twin_lj_fswitch_noener, nbnxn_kernel_ref_tab_twin_lj_pswitch_noener }
};

p_nbk_func_ener p_nbk_c_ener[coultNR][vdwtNR] =
{
    { nbnxn_kernel_ref_rf_lj_cut_ener,       nbnxn_kernel_ref_rf_lj_fswitch_ener,       nbnxn_kernel_ref_rf_lj_pswitch_ener },
    { nbnxn_kernel_ref_tab_lj_cut_ener,      nbnxn_kernel_ref_tab_lj_fswitch_ener,      nbnxn_kernel_ref_tab_lj_pswitch_ener },
    { nbnxn_kernel_ref_tab_twin_lj_cut_ener, nbnxn_kernel_ref_tab_twin_lj_fswitch_ener, nbnxn_kernel_ref_tab_twin_lj_pswitch_ener }
};

p_nbk_func_ener p_nbk_c_energrp[coultNR][vdwtNR] =
{
    { nbnxn_kernel_ref_rf_lj_cut_energrp,       nbnxn_kernel_ref_rf_lj_fswitch_energrp,       nbnxn_kernel_ref_rf_lj_pswitch_energrp },
    { nbnxn_kernel_ref_tab_lj_cut_energrp,      nbnxn_kernel_ref_tab_lj_fswitch_energrp,      nbnxn_kernel_ref_tab_lj_pswitch_energrp },
    { nbnxn_kernel_ref_tab_twin_lj_cut_energrp, nbnxn_kernel_ref_tab_twin_lj_fswitch_energrp, nbnxn_kernel_ref_tab_twin_lj_pswitch_energrp }
};

void
nbnxn_kernel_ref(const nbnxn_pairlist_set_t *nbl_list,
                 const nbnxn_atomdata_t     *nbat,
                 const interaction_const_t  *ic,
                 rvec                       *shift_vec,
                 int                         force_flags,
                 int                         clearF,
                 real                       *fshift,
                 real                       *Vc,
                 real                       *Vvdw)
{
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                coult;
    int                vdwt;
    int                nb;

    nnbl = nbl_list->nnbl;
    nbl  = nbl_list->nbl;

    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        if (ic->rcoulomb == ic->rvdw)
        {
            coult = coultTAB;
        }
        else
        {
            coult = coultTAB_TWIN;
        }
    }

    switch (ic->vdw_modifier)
    {
        case eintmodPOTSHIFT:
        case eintmodNONE:
            vdwt = vdwtCUT;
            break;
        case eintmodFORCESWITCH:
            vdwt = vdwtFSWITCH;
            break;
        case eintmodPOTSWITCH:
            vdwt = vdwtPSWITCH;
            break;
        default:
            gmx_incons("Unsupported VdW modifier");
            break;
    }

#pragma omp parallel for schedule(static) num_threads(gmx_omp_nthreads_get(emntNonbonded))
    for (nb = 0; nb < nnbl; nb++)
    {
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
            p_nbk_c_noener[coult][vdwt](nbl[nb], nbat,
                                        ic,
                                        shift_vec,
                                        out->f,
                                        fshift_p);
        }
        else if (out->nV == 1)
        {
            /* No energy groups */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            p_nbk_c_ener[coult][vdwt](nbl[nb], nbat,
                                      ic,
                                      shift_vec,
                                      out->f,
                                      fshift_p,
                                      out->Vvdw,
                                      out->Vc);
        }
        else
        {
            /* Calculate energy group contributions */
            int i;

            for (i = 0; i < out->nV; i++)
            {
                out->Vvdw[i] = 0;
            }
            for (i = 0; i < out->nV; i++)
            {
                out->Vc[i] = 0;
            }

            p_nbk_c_energrp[coult][vdwt](nbl[nb], nbat,
                                         ic,
                                         shift_vec,
                                         out->f,
                                         fshift_p,
                                         out->Vvdw,
                                         out->Vc);
        }
    }

    if (force_flags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(nbat, nnbl, Vvdw, Vc);
    }
}
