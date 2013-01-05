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

/* GMX_MM128_HERE or GMX_MM256_HERE should be set before including this file.
 * gmx_sse_or_avh.h should be included before including this file.
 */

/* Copies PBC shifted i-cell packed atom coordinates to working array */
#ifdef GMX_MM128_HERE
static void icell_set_x_x86_simd128
#else
#ifdef GMX_MM256_HERE
static void icell_set_x_x86_simd256
#else
"error: GMX_MM128_HERE or GMX_MM256_HERE not defined"
#endif
#endif
(int ci,
 real shx, real shy, real shz,
 int na_c,
 int stride, const real *x,
 nbnxn_list_work_t *work)
{
    int  ia;
#ifdef GMX_MM128_HERE
    nbnxn_x_ci_x86_simd128_t *x_ci;

    x_ci = work->x_ci_x86_simd128;

    ia = X_IND_CI_S128(ci);
#else
    nbnxn_x_ci_x86_simd256_t *x_ci;

    x_ci = work->x_ci_x86_simd256;

    ia = X_IND_CI_S256(ci);
#endif

    x_ci->ix_SSE0 = gmx_set1_pr(x[ia + 0*STRIDE_S    ] + shx);
    x_ci->iy_SSE0 = gmx_set1_pr(x[ia + 1*STRIDE_S    ] + shy);
    x_ci->iz_SSE0 = gmx_set1_pr(x[ia + 2*STRIDE_S    ] + shz);
    x_ci->ix_SSE1 = gmx_set1_pr(x[ia + 0*STRIDE_S + 1] + shx);
    x_ci->iy_SSE1 = gmx_set1_pr(x[ia + 1*STRIDE_S + 1] + shy);
    x_ci->iz_SSE1 = gmx_set1_pr(x[ia + 2*STRIDE_S + 1] + shz);
    x_ci->ix_SSE2 = gmx_set1_pr(x[ia + 0*STRIDE_S + 2] + shx);
    x_ci->iy_SSE2 = gmx_set1_pr(x[ia + 1*STRIDE_S + 2] + shy);
    x_ci->iz_SSE2 = gmx_set1_pr(x[ia + 2*STRIDE_S + 2] + shz);
    x_ci->ix_SSE3 = gmx_set1_pr(x[ia + 0*STRIDE_S + 3] + shx);
    x_ci->iy_SSE3 = gmx_set1_pr(x[ia + 1*STRIDE_S + 3] + shy);
    x_ci->iz_SSE3 = gmx_set1_pr(x[ia + 2*STRIDE_S + 3] + shz);
}

/* SSE or AVX code for making a pair list of cell ci vs cell cjf-cjl
 * for coordinates in packed format.
 * Checks bouding box distances and possibly atom pair distances.
 * This is an accelerated version of make_cluster_list_simple.
 */
#ifdef GMX_MM128_HERE
static void make_cluster_list_x86_simd128
#else
#ifdef GMX_MM256_HERE
static void make_cluster_list_x86_simd256
#else
"error: GMX_MM128_HERE or GMX_MM256_HERE not defined"
#endif
#endif
(const nbnxn_grid_t *gridj,
 nbnxn_pairlist_t *nbl,
 int ci, int cjf, int cjl,
 gmx_bool remove_sub_diag,
 const real *x_j,
 real rl2, float rbb2,
 int *ndistc)
{
#ifdef GMX_MM128_HERE
    const nbnxn_x_ci_x86_simd128_t *work;
#else
    const nbnxn_x_ci_x86_simd256_t *work;
#endif

    const float *bb_ci;

    gmx_mm_pr    jx_SSE, jy_SSE, jz_SSE;

    gmx_mm_pr    dx_SSE0, dy_SSE0, dz_SSE0;
    gmx_mm_pr    dx_SSE1, dy_SSE1, dz_SSE1;
    gmx_mm_pr    dx_SSE2, dy_SSE2, dz_SSE2;
    gmx_mm_pr    dx_SSE3, dy_SSE3, dz_SSE3;

    gmx_mm_pr    rsq_SSE0;
    gmx_mm_pr    rsq_SSE1;
    gmx_mm_pr    rsq_SSE2;
    gmx_mm_pr    rsq_SSE3;

    gmx_mm_pr    wco_SSE0;
    gmx_mm_pr    wco_SSE1;
    gmx_mm_pr    wco_SSE2;
    gmx_mm_pr    wco_SSE3;
    gmx_mm_pr    wco_any_SSE01, wco_any_SSE23, wco_any_SSE;

    gmx_mm_pr    rc2_SSE;

    gmx_bool     InRange;
    float        d2;
    int          xind_f, xind_l, cj;

#ifdef GMX_MM128_HERE
    cjf = CI_TO_CJ_S128(cjf);
    cjl = CI_TO_CJ_S128(cjl+1) - 1;

    work = nbl->work->x_ci_x86_simd128;
#else
    cjf = CI_TO_CJ_S256(cjf);
    cjl = CI_TO_CJ_S256(cjl+1) - 1;

    work = nbl->work->x_ci_x86_simd256;
#endif

    bb_ci = nbl->work->bb_ci;

    rc2_SSE   = gmx_set1_pr(rl2);

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
        d2       = subc_bb_dist2_sse(4, 0, bb_ci, cjf, gridj->bbj);
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
#ifdef GMX_MM128_HERE
            xind_f  = X_IND_CJ_S128(CI_TO_CJ_S128(gridj->cell0) + cjf);
#else
            xind_f  = X_IND_CJ_S256(CI_TO_CJ_S256(gridj->cell0) + cjf);
#endif
            jx_SSE  = gmx_load_pr(x_j+xind_f+0*STRIDE_S);
            jy_SSE  = gmx_load_pr(x_j+xind_f+1*STRIDE_S);
            jz_SSE  = gmx_load_pr(x_j+xind_f+2*STRIDE_S);


            /* Calculate distance */
            dx_SSE0            = gmx_sub_pr(work->ix_SSE0, jx_SSE);
            dy_SSE0            = gmx_sub_pr(work->iy_SSE0, jy_SSE);
            dz_SSE0            = gmx_sub_pr(work->iz_SSE0, jz_SSE);
            dx_SSE1            = gmx_sub_pr(work->ix_SSE1, jx_SSE);
            dy_SSE1            = gmx_sub_pr(work->iy_SSE1, jy_SSE);
            dz_SSE1            = gmx_sub_pr(work->iz_SSE1, jz_SSE);
            dx_SSE2            = gmx_sub_pr(work->ix_SSE2, jx_SSE);
            dy_SSE2            = gmx_sub_pr(work->iy_SSE2, jy_SSE);
            dz_SSE2            = gmx_sub_pr(work->iz_SSE2, jz_SSE);
            dx_SSE3            = gmx_sub_pr(work->ix_SSE3, jx_SSE);
            dy_SSE3            = gmx_sub_pr(work->iy_SSE3, jy_SSE);
            dz_SSE3            = gmx_sub_pr(work->iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_calc_rsq_pr(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_calc_rsq_pr(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_calc_rsq_pr(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_calc_rsq_pr(dx_SSE3, dy_SSE3, dz_SSE3);

            wco_SSE0           = gmx_cmplt_pr(rsq_SSE0, rc2_SSE);
            wco_SSE1           = gmx_cmplt_pr(rsq_SSE1, rc2_SSE);
            wco_SSE2           = gmx_cmplt_pr(rsq_SSE2, rc2_SSE);
            wco_SSE3           = gmx_cmplt_pr(rsq_SSE3, rc2_SSE);

            wco_any_SSE01      = gmx_or_pr(wco_SSE0, wco_SSE1);
            wco_any_SSE23      = gmx_or_pr(wco_SSE2, wco_SSE3);
            wco_any_SSE        = gmx_or_pr(wco_any_SSE01, wco_any_SSE23);

            InRange            = gmx_movemask_pr(wco_any_SSE);

            *ndistc += 4*GMX_X86_SIMD_WIDTH_HERE;
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
        d2       = subc_bb_dist2_sse(4, 0, bb_ci, cjl, gridj->bbj);
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
#ifdef GMX_MM128_HERE
            xind_l  = X_IND_CJ_S128(CI_TO_CJ_S128(gridj->cell0) + cjl);
#else
            xind_l  = X_IND_CJ_S256(CI_TO_CJ_S256(gridj->cell0) + cjl);
#endif
            jx_SSE  = gmx_load_pr(x_j+xind_l+0*STRIDE_S);
            jy_SSE  = gmx_load_pr(x_j+xind_l+1*STRIDE_S);
            jz_SSE  = gmx_load_pr(x_j+xind_l+2*STRIDE_S);

            /* Calculate distance */
            dx_SSE0            = gmx_sub_pr(work->ix_SSE0, jx_SSE);
            dy_SSE0            = gmx_sub_pr(work->iy_SSE0, jy_SSE);
            dz_SSE0            = gmx_sub_pr(work->iz_SSE0, jz_SSE);
            dx_SSE1            = gmx_sub_pr(work->ix_SSE1, jx_SSE);
            dy_SSE1            = gmx_sub_pr(work->iy_SSE1, jy_SSE);
            dz_SSE1            = gmx_sub_pr(work->iz_SSE1, jz_SSE);
            dx_SSE2            = gmx_sub_pr(work->ix_SSE2, jx_SSE);
            dy_SSE2            = gmx_sub_pr(work->iy_SSE2, jy_SSE);
            dz_SSE2            = gmx_sub_pr(work->iz_SSE2, jz_SSE);
            dx_SSE3            = gmx_sub_pr(work->ix_SSE3, jx_SSE);
            dy_SSE3            = gmx_sub_pr(work->iy_SSE3, jy_SSE);
            dz_SSE3            = gmx_sub_pr(work->iz_SSE3, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_calc_rsq_pr(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE1           = gmx_calc_rsq_pr(dx_SSE1, dy_SSE1, dz_SSE1);
            rsq_SSE2           = gmx_calc_rsq_pr(dx_SSE2, dy_SSE2, dz_SSE2);
            rsq_SSE3           = gmx_calc_rsq_pr(dx_SSE3, dy_SSE3, dz_SSE3);

            wco_SSE0           = gmx_cmplt_pr(rsq_SSE0, rc2_SSE);
            wco_SSE1           = gmx_cmplt_pr(rsq_SSE1, rc2_SSE);
            wco_SSE2           = gmx_cmplt_pr(rsq_SSE2, rc2_SSE);
            wco_SSE3           = gmx_cmplt_pr(rsq_SSE3, rc2_SSE);

            wco_any_SSE01      = gmx_or_pr(wco_SSE0, wco_SSE1);
            wco_any_SSE23      = gmx_or_pr(wco_SSE2, wco_SSE3);
            wco_any_SSE        = gmx_or_pr(wco_any_SSE01, wco_any_SSE23);

            InRange            = gmx_movemask_pr(wco_any_SSE);

            *ndistc += 4*GMX_X86_SIMD_WIDTH_HERE;
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
#ifdef GMX_MM128_HERE
            nbl->cj[nbl->ncj].cj   = CI_TO_CJ_S128(gridj->cell0) + cj;
            nbl->cj[nbl->ncj].excl = get_imask_x86_simd128(remove_sub_diag, ci, cj);
#else
            nbl->cj[nbl->ncj].cj   = CI_TO_CJ_S256(gridj->cell0) + cj;
            nbl->cj[nbl->ncj].excl = get_imask_x86_simd256(remove_sub_diag, ci, cj);
#endif
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}
