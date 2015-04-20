/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

/* Get the half-width SIMD stuff from the kernel utils files */
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils.h"


#if GMX_SIMD_REAL_WIDTH >= 2*NBNXN_CPU_CLUSTER_I_SIZE
#define STRIDE_S  (GMX_SIMD_REAL_WIDTH/2)
#else
#define STRIDE_S  NBNXN_CPU_CLUSTER_I_SIZE
#endif

static gmx_inline gmx_simd_real_t gmx_load_hpr_hilo_pr(const real *a)
{
    gmx_mm_hpr       a_S;
    gmx_simd_real_t  a_a_S;

    gmx_load_hpr(&a_S, a);

    gmx_2hpr_to_pr(a_S, a_S, &a_a_S);

    return a_a_S;
}

static gmx_inline gmx_simd_real_t gmx_set_2real_shift_pr(const real *a, real shift)
{
    gmx_mm_hpr       a0_S, a1_S;
    gmx_simd_real_t  a0_a1_S;

    gmx_set1_hpr(&a0_S, a[0] + shift);
    gmx_set1_hpr(&a1_S, a[1] + shift);

    gmx_2hpr_to_pr(a0_S, a1_S, &a0_a1_S);

    return a0_a1_S;
}

/* Copies PBC shifted i-cell packed atom coordinates to working array */
static gmx_inline void
icell_set_x_simd_2xnn(int ci,
                      real shx, real shy, real shz,
                      int gmx_unused na_c,
                      int gmx_unused stride, const real *x,
                      nbnxn_list_work_t *work)
{
    int                     ia;
    nbnxn_x_ci_simd_2xnn_t *x_ci;

    x_ci = work->x_ci_simd_2xnn;

    ia = X_IND_CI_SIMD_2XNN(ci);

    x_ci->ix_S0 = gmx_set_2real_shift_pr(x + ia + 0*STRIDE_S + 0, shx);
    x_ci->iy_S0 = gmx_set_2real_shift_pr(x + ia + 1*STRIDE_S + 0, shy);
    x_ci->iz_S0 = gmx_set_2real_shift_pr(x + ia + 2*STRIDE_S + 0, shz);
    x_ci->ix_S2 = gmx_set_2real_shift_pr(x + ia + 0*STRIDE_S + 2, shx);
    x_ci->iy_S2 = gmx_set_2real_shift_pr(x + ia + 1*STRIDE_S + 2, shy);
    x_ci->iz_S2 = gmx_set_2real_shift_pr(x + ia + 2*STRIDE_S + 2, shz);
}

/* SIMD code for making a pair list of cell ci vs cell cjf-cjl
 * for coordinates in packed format.
 * Checks bouding box distances and possibly atom pair distances.
 * This is an accelerated version of make_cluster_list_simple.
 */
static gmx_inline void
make_cluster_list_simd_2xnn(const nbnxn_grid_t *gridj,
                            nbnxn_pairlist_t *nbl,
                            int ci, int cjf, int cjl,
                            gmx_bool remove_sub_diag,
                            const real *x_j,
                            real rl2, float rbb2,
                            int *ndistc)
{
    const nbnxn_x_ci_simd_2xnn_t       *work;
    const nbnxn_bb_t                   *bb_ci;

    gmx_simd_real_t                     jx_S, jy_S, jz_S;

    gmx_simd_real_t                     dx_S0, dy_S0, dz_S0;
    gmx_simd_real_t                     dx_S2, dy_S2, dz_S2;

    gmx_simd_real_t                     rsq_S0;
    gmx_simd_real_t                     rsq_S2;

    gmx_simd_bool_t                     wco_S0;
    gmx_simd_bool_t                     wco_S2;
    gmx_simd_bool_t                     wco_any_S;

    gmx_simd_real_t                     rc2_S;

    gmx_bool                            InRange;
    float                               d2;
    int                                 xind_f, xind_l, cj;

    cjf = CI_TO_CJ_SIMD_2XNN(cjf);
    cjl = CI_TO_CJ_SIMD_2XNN(cjl+1) - 1;

    work = nbl->work->x_ci_simd_2xnn;

    bb_ci = nbl->work->bb_ci;

    rc2_S   = gmx_simd_set1_r(rl2);

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
#ifdef NBNXN_SEARCH_BB_SIMD4
        d2 = subc_bb_dist2_simd4(0, bb_ci, cjf, gridj->bbj);
#else
        d2 = subc_bb_dist2(0, bb_ci, cjf, gridj->bbj);
#endif
        *ndistc += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rl2)
        {
            xind_f  = X_IND_CJ_SIMD_2XNN(CI_TO_CJ_SIMD_2XNN(gridj->cell0) + cjf);

            jx_S  = gmx_load_hpr_hilo_pr(x_j+xind_f+0*STRIDE_S);
            jy_S  = gmx_load_hpr_hilo_pr(x_j+xind_f+1*STRIDE_S);
            jz_S  = gmx_load_hpr_hilo_pr(x_j+xind_f+2*STRIDE_S);

            /* Calculate distance */
            dx_S0            = gmx_simd_sub_r(work->ix_S0, jx_S);
            dy_S0            = gmx_simd_sub_r(work->iy_S0, jy_S);
            dz_S0            = gmx_simd_sub_r(work->iz_S0, jz_S);
            dx_S2            = gmx_simd_sub_r(work->ix_S2, jx_S);
            dy_S2            = gmx_simd_sub_r(work->iy_S2, jy_S);
            dz_S2            = gmx_simd_sub_r(work->iz_S2, jz_S);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_S0           = gmx_simd_calc_rsq_r(dx_S0, dy_S0, dz_S0);
            rsq_S2           = gmx_simd_calc_rsq_r(dx_S2, dy_S2, dz_S2);

            wco_S0           = gmx_simd_cmplt_r(rsq_S0, rc2_S);
            wco_S2           = gmx_simd_cmplt_r(rsq_S2, rc2_S);

            wco_any_S        = gmx_simd_or_b(wco_S0, wco_S2);

            InRange          = gmx_simd_anytrue_b(wco_any_S);

            *ndistc += 2*GMX_SIMD_REAL_WIDTH;
        }
        if (!InRange)
        {
            cjf++;
        }
    }
    if (!InRange)
    {
        return;
    }

    InRange = FALSE;
    while (!InRange && cjl > cjf)
    {
#ifdef NBNXN_SEARCH_BB_SIMD4
        d2 = subc_bb_dist2_simd4(0, bb_ci, cjl, gridj->bbj);
#else
        d2 = subc_bb_dist2(0, bb_ci, cjl, gridj->bbj);
#endif
        *ndistc += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rl2)
        {
            xind_l  = X_IND_CJ_SIMD_2XNN(CI_TO_CJ_SIMD_2XNN(gridj->cell0) + cjl);

            jx_S  = gmx_load_hpr_hilo_pr(x_j+xind_l+0*STRIDE_S);
            jy_S  = gmx_load_hpr_hilo_pr(x_j+xind_l+1*STRIDE_S);
            jz_S  = gmx_load_hpr_hilo_pr(x_j+xind_l+2*STRIDE_S);

            /* Calculate distance */
            dx_S0            = gmx_simd_sub_r(work->ix_S0, jx_S);
            dy_S0            = gmx_simd_sub_r(work->iy_S0, jy_S);
            dz_S0            = gmx_simd_sub_r(work->iz_S0, jz_S);
            dx_S2            = gmx_simd_sub_r(work->ix_S2, jx_S);
            dy_S2            = gmx_simd_sub_r(work->iy_S2, jy_S);
            dz_S2            = gmx_simd_sub_r(work->iz_S2, jz_S);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_S0           = gmx_simd_calc_rsq_r(dx_S0, dy_S0, dz_S0);
            rsq_S2           = gmx_simd_calc_rsq_r(dx_S2, dy_S2, dz_S2);

            wco_S0           = gmx_simd_cmplt_r(rsq_S0, rc2_S);
            wco_S2           = gmx_simd_cmplt_r(rsq_S2, rc2_S);

            wco_any_S        = gmx_simd_or_b(wco_S0, wco_S2);

            InRange          = gmx_simd_anytrue_b(wco_any_S);

            *ndistc += 2*GMX_SIMD_REAL_WIDTH;
        }
        if (!InRange)
        {
            cjl--;
        }
    }

    if (cjf <= cjl)
    {
        for (cj = cjf; cj <= cjl; cj++)
        {
            /* Store cj and the interaction mask */
            nbl->cj[nbl->ncj].cj   = CI_TO_CJ_SIMD_2XNN(gridj->cell0) + cj;
            nbl->cj[nbl->ncj].excl = get_imask_simd_2xnn(remove_sub_diag, ci, cj);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

#undef STRIDE_S
