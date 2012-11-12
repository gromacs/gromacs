/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
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
#include "../nbnxn_consts.h"
#include "nbnxn_kernel_common.h"

#ifdef GMX_X86_SSE2

#include "nbnxn_kernel_x86_simd128.h"

/* Include all flavors of the 128-bit SSE or AVX kernel loops */

#define GMX_MM128_HERE

/* Analytical reaction-field kernels */
#define CALC_COUL_RF

#include "nbnxn_kernel_x86_simd_includes.h"

#undef CALC_COUL_RF

/* Tabulated exclusion interaction electrostatics kernels */
#define CALC_COUL_TAB

/* Single cut-off: rcoulomb = rvdw */
#include "nbnxn_kernel_x86_simd_includes.h"

/* Twin cut-off: rcoulomb >= rvdw */
#define VDW_CUTOFF_CHECK
#include "nbnxn_kernel_x86_simd_includes.h"
#undef VDW_CUTOFF_CHECK

#undef CALC_COUL_TAB


typedef void (*p_nbk_func_ener)(const nbnxn_pairlist_t *   nbl,
                                const nbnxn_atomdata_t *   nbat,
                                const interaction_const_t *ic,
                                rvec *                     shift_vec,
                                real *                     f,
                                real *                     fshift,
                                real *                     Vvdw,
                                real *                     Vc);

typedef void (*p_nbk_func_noener)(const nbnxn_pairlist_t *   nbl,
                                  const nbnxn_atomdata_t *   nbat,
                                  const interaction_const_t *ic,
                                  rvec *                     shift_vec,
                                  real *                     f,
                                  real *                     fshift);

enum {
    coultRF, coultTAB, coultTAB_TWIN, coultNR
};


static p_nbk_func_ener p_nbk_ener[coultNR][ljcrNR] =
{ { nbnxn_kernel_x86_simd128_rf_comb_geom_ener,
    nbnxn_kernel_x86_simd128_rf_comb_lb_ener,
    nbnxn_kernel_x86_simd128_rf_comb_none_ener },
  { nbnxn_kernel_x86_simd128_tab_comb_geom_ener,
    nbnxn_kernel_x86_simd128_tab_comb_lb_ener,
    nbnxn_kernel_x86_simd128_tab_twin_comb_none_ener },
  { nbnxn_kernel_x86_simd128_tab_twin_comb_geom_ener,
    nbnxn_kernel_x86_simd128_tab_twin_comb_lb_ener,
    nbnxn_kernel_x86_simd128_tab_twin_comb_none_ener }  };

static p_nbk_func_ener p_nbk_energrp[coultNR][ljcrNR] =
{ { nbnxn_kernel_x86_simd128_rf_comb_geom_energrp,
    nbnxn_kernel_x86_simd128_rf_comb_lb_energrp,
    nbnxn_kernel_x86_simd128_rf_comb_none_energrp },
  { nbnxn_kernel_x86_simd128_tab_comb_geom_energrp,
    nbnxn_kernel_x86_simd128_tab_comb_lb_energrp,
    nbnxn_kernel_x86_simd128_tab_comb_none_energrp },
  { nbnxn_kernel_x86_simd128_tab_twin_comb_geom_energrp,
    nbnxn_kernel_x86_simd128_tab_twin_comb_lb_energrp,
    nbnxn_kernel_x86_simd128_tab_twin_comb_none_energrp } };

static p_nbk_func_noener p_nbk_noener[coultNR][ljcrNR] =
{ { nbnxn_kernel_x86_simd128_rf_comb_geom_noener,
    nbnxn_kernel_x86_simd128_rf_comb_lb_noener,
    nbnxn_kernel_x86_simd128_rf_comb_none_noener },
  { nbnxn_kernel_x86_simd128_tab_comb_geom_noener,
    nbnxn_kernel_x86_simd128_tab_comb_lb_noener,
    nbnxn_kernel_x86_simd128_tab_comb_none_noener },
  { nbnxn_kernel_x86_simd128_tab_twin_comb_geom_noener,
    nbnxn_kernel_x86_simd128_tab_twin_comb_lb_noener,
    nbnxn_kernel_x86_simd128_tab_twin_comb_none_noener } };


static void reduce_group_energies(int ng, int ng_2log,
                                  const real *VSvdw, const real *VSc,
                                  real *Vvdw, real *Vc)
{
    int ng_p2, i, j, j0, j1, c, s;

#define SIMD_WIDTH       (GMX_X86_SIMD_WIDTH_HERE)
#define SIMD_WIDTH_HALF  (GMX_X86_SIMD_WIDTH_HERE / 2)

    ng_p2 = (1 << ng_2log);

    /* The size of the x86 SIMD energy group buffer array is:
     * ng*ng*ng_p2*SIMD_WIDTH_HALF*SIMD_WIDTH
     */
    for(i = 0; i < ng; i++)
    {
        for(j = 0; j < ng; j++)
        {
            Vvdw[i * ng + j] = 0;
            Vc[i * ng + j]   = 0;
        }

        for(j1 = 0; j1 < ng; j1++)
        {
            for(j0 = 0; j0 < ng; j0++)
            {
                c = ((i * ng + j1) * ng_p2 + j0) * SIMD_WIDTH_HALF * SIMD_WIDTH;
                for(s = 0; s < SIMD_WIDTH_HALF; s++)
                {
                    Vvdw[i * ng + j0] += VSvdw[c + 0];
                    Vvdw[i * ng + j1] += VSvdw[c + 1];
                    Vc  [i * ng + j0] += VSc  [c + 0];
                    Vc  [i * ng + j1] += VSc  [c + 1];
                    c                 += SIMD_WIDTH + 2;
                }
            }
        }
    }
}

#endif /* GMX_X86_SSE2 */

void
nbnxn_kernel_x86_simd128(nbnxn_pairlist_set_t *     nbl_list,
                         const nbnxn_atomdata_t *   nbat,
                         const interaction_const_t *ic,
                         rvec *                     shift_vec,
                         int                        force_flags,
                         int                        clearF,
                         real *                     fshift,
                         real *                     Vc,
                         real *                     Vvdw)
#ifdef GMX_X86_SSE2
{
    int nnbl;
    nbnxn_pairlist_t **nbl;
    int coult;
    int nb;

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

#pragma omp parallel for schedule(static) num_threads(gmx_omp_nthreads_get(emntNonbonded))
    for(nb = 0; nb < nnbl; nb++)
    {
        nbnxn_atomdata_output_t *out;
        real *fshift_p;

        out = &nbat->out[nb];

        if (clearF == enbvClearFYes)
        {
            clear_f(nbat, out->f);
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

        /* With Ewald type electrostatics we the forces for excluded atom pairs
         * should not contribute to the virial sum. The exclusion forces
         * are not calculate in the energy kernels, but are in _noener.
         */
        if (!((force_flags & GMX_FORCE_ENERGY) ||
              (EEL_FULL(ic->eeltype) && (force_flags & GMX_FORCE_VIRIAL))))
        {
            /* Don't calculate energies */
            p_nbk_noener[coult][nbat->comb_rule](nbl[nb], nbat,
                                                 ic,
                                                 shift_vec,
                                                 out->f,
                                                 fshift_p);
        }
        else if (out->nV == 1 || !(force_flags & GMX_FORCE_ENERGY))
        {
            /* No energy groups */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            p_nbk_ener[coult][nbat->comb_rule](nbl[nb], nbat,
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

            for(i = 0; i < out->nVS; i++)
            {
                out->VSvdw[i] = 0;
            }
            for(i = 0; i < out->nVS; i++)
            {
                out->VSc[i] = 0;
            }

            p_nbk_energrp[coult][nbat->comb_rule](nbl[nb], nbat,
                                                  ic,
                                                  shift_vec,
                                                  out->f,
                                                  fshift_p,
                                                  out->VSvdw,
                                                  out->VSc);

            reduce_group_energies(nbat->nenergrp, nbat->neg_2log,
                                  out->VSvdw, out->VSc,
                                  out->Vvdw, out->Vc);
        }
    }

    if (force_flags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(nbat, nnbl, Vvdw, Vc);
    }
}
#else
{
    gmx_incons("nbnxn_kernel_x86_simd128 called while GROMACS was configured without SSE enabled");
}
#endif
