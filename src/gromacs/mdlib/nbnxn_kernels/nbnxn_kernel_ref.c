/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
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
