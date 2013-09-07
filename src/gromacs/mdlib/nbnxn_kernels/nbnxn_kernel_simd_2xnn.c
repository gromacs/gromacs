/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include "../nbnxn_consts.h"
#include "nbnxn_kernel_common.h"

#ifdef GMX_NBNXN_SIMD_2XNN

/* Include the full width SIMD macros */
#include "gmx_simd_macros.h"
#include "gmx_simd_vec.h"

#include "nbnxn_kernel_simd_2xnn.h"

#if !(GMX_SIMD_WIDTH_HERE == 8 || GMX_SIMD_WIDTH_HERE == 16)
#error "unsupported SIMD width"
#endif

#define SUM_SIMD4(x) (x[0]+x[1]+x[2]+x[3])

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    (GMX_SIMD_WIDTH_HERE/2)

/* The stride of all the atom data arrays is equal to half the SIMD width */
#define STRIDE     (GMX_SIMD_WIDTH_HERE/2)

#if GMX_SIMD_WIDTH_HERE == 8
#define SUM_SIMD(x) (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7])
#else
#if GMX_SIMD_WIDTH_HERE == 16
/* This is getting ridiculous, SIMD horizontal adds would help,
 * but this is not performance critical (only used to reduce energies)
 */
#define SUM_SIMD(x) (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15])
#else
#error "unsupported kernel configuration"
#endif
#endif


#include "nbnxn_kernel_simd_utils.h"

static inline void
gmx_load_simd_2xnn_interactions(int            excl,
                                gmx_exclfilter filter_S0,
                                gmx_exclfilter filter_S2,
                                gmx_mm_pb     *interact_S0,
                                gmx_mm_pb     *interact_S2)
{
    /* Load integer topology exclusion interaction mask */
    gmx_exclfilter mask_pr_S = gmx_load1_exclfilter(excl);
    *interact_S0  = gmx_checkbitmask_pb(mask_pr_S, filter_S0);
    *interact_S2  = gmx_checkbitmask_pb(mask_pr_S, filter_S2);
}

/* Include all flavors of the SSE or AVX 2x(N+N) kernel loops */

/* Analytical reaction-field kernels */
#define CALC_COUL_RF

#include "nbnxn_kernel_simd_2xnn_includes.h"

#undef CALC_COUL_RF

/* Tabulated exclusion interaction electrostatics kernels */
#define CALC_COUL_TAB

/* Single cut-off: rcoulomb = rvdw */
#include "nbnxn_kernel_simd_2xnn_includes.h"

/* Twin cut-off: rcoulomb >= rvdw */
#define VDW_CUTOFF_CHECK
#include "nbnxn_kernel_simd_2xnn_includes.h"
#undef VDW_CUTOFF_CHECK

#undef CALC_COUL_TAB

/* Analytical Ewald exclusion interaction electrostatics kernels */
#define CALC_COUL_EWALD

/* Single cut-off: rcoulomb = rvdw */
#include "nbnxn_kernel_simd_2xnn_includes.h"

/* Twin cut-off: rcoulomb >= rvdw */
#define VDW_CUTOFF_CHECK
#include "nbnxn_kernel_simd_2xnn_includes.h"
#undef VDW_CUTOFF_CHECK

#undef CALC_COUL_EWALD


typedef void (*p_nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                                const nbnxn_atomdata_t     *nbat,
                                const interaction_const_t  *ic,
                                rvec                       *shift_vec,
                                real                       *f,
                                real                       *fshift,
                                real                       *Vvdw,
                                real                       *Vc);

typedef void (*p_nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                                  const nbnxn_atomdata_t     *nbat,
                                  const interaction_const_t  *ic,
                                  rvec                       *shift_vec,
                                  real                       *f,
                                  real                       *fshift);

enum {
    coultRF, coultTAB, coultTAB_TWIN, coultEWALD, coultEWALD_TWIN, coultNR
};

#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _ener
static p_nbk_func_ener p_nbk_ener[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom), NBK_FN(rf, lb), NBK_FN(rf, none) },
  { NBK_FN(tab, geom), NBK_FN(tab, lb), NBK_FN(tab, none) },
  { NBK_FN(tab_twin, geom), NBK_FN(tab_twin, lb), NBK_FN(tab_twin, none) },
  { NBK_FN(ewald, geom), NBK_FN(ewald, lb), NBK_FN(ewald, none) },
  { NBK_FN(ewald_twin, geom), NBK_FN(ewald_twin, lb), NBK_FN(ewald_twin, none) } };
#undef NBK_FN

#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _energrp
static p_nbk_func_ener p_nbk_energrp[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom), NBK_FN(rf, lb), NBK_FN(rf, none) },
  { NBK_FN(tab, geom), NBK_FN(tab, lb), NBK_FN(tab, none) },
  { NBK_FN(tab_twin, geom), NBK_FN(tab_twin, lb), NBK_FN(tab_twin, none) },
  { NBK_FN(ewald, geom), NBK_FN(ewald, lb), NBK_FN(ewald, none) },
  { NBK_FN(ewald_twin, geom), NBK_FN(ewald_twin, lb), NBK_FN(ewald_twin, none) } };
#undef NBK_FN

#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _noener
static p_nbk_func_noener p_nbk_noener[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom), NBK_FN(rf, lb), NBK_FN(rf, none) },
  { NBK_FN(tab, geom), NBK_FN(tab, lb), NBK_FN(tab, none) },
  { NBK_FN(tab_twin, geom), NBK_FN(tab_twin, lb), NBK_FN(tab_twin, none) },
  { NBK_FN(ewald, geom), NBK_FN(ewald, lb), NBK_FN(ewald, none) },
  { NBK_FN(ewald_twin, geom), NBK_FN(ewald_twin, lb), NBK_FN(ewald_twin, none) } };
#undef NBK_FN


static void reduce_group_energies(int ng, int ng_2log,
                                  const real *VSvdw, const real *VSc,
                                  real *Vvdw, real *Vc)
{
    const int unrollj      = GMX_SIMD_WIDTH_HERE/2;
    const int unrollj_half = unrollj/2;
    int       ng_p2, i, j, j0, j1, c, s;

    ng_p2 = (1<<ng_2log);

    /* The size of the x86 SIMD energy group buffer array is:
     * ng*ng*ng_p2*unrollj_half*simd_width
     */
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            Vvdw[i*ng+j] = 0;
            Vc[i*ng+j]   = 0;
        }

        for (j1 = 0; j1 < ng; j1++)
        {
            for (j0 = 0; j0 < ng; j0++)
            {
                c = ((i*ng + j1)*ng_p2 + j0)*unrollj_half*unrollj;
                for (s = 0; s < unrollj_half; s++)
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

#endif /* GMX_NBNXN_SIMD_2XNN */

void
nbnxn_kernel_simd_2xnn(nbnxn_pairlist_set_t       gmx_unused *nbl_list,
                       const nbnxn_atomdata_t     gmx_unused *nbat,
                       const interaction_const_t  gmx_unused *ic,
                       int                        gmx_unused ewald_excl,
                       rvec                       gmx_unused *shift_vec,
                       int                        gmx_unused  force_flags,
                       int                        gmx_unused  clearF,
                       real                       gmx_unused *fshift,
                       real                       gmx_unused *Vc,
                       real                       gmx_unused *Vvdw)
#ifdef GMX_NBNXN_SIMD_2XNN
{
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                coult;
    int                nb;

    nnbl = nbl_list->nnbl;
    nbl  = nbl_list->nbl;

    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        if (ewald_excl == ewaldexclTable)
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
        else
        {
            if (ic->rcoulomb == ic->rvdw)
            {
                coult = coultEWALD;
            }
            else
            {
                coult = coultEWALD_TWIN;
            }
        }
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

            for (i = 0; i < out->nVS; i++)
            {
                out->VSvdw[i] = 0;
            }
            for (i = 0; i < out->nVS; i++)
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
    gmx_incons("nbnxn_kernel_simd_2xnn called while GROMACS was configured without 2x(N+N) SIMD kernels enabled");
}
#endif
