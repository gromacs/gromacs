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

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

#if GMX_SIMD_REAL_WIDTH >= NBNXN_CPU_CLUSTER_I_SIZE
#define STRIDE_S  (GMX_SIMD_REAL_WIDTH)
#else
#define STRIDE_S  NBNXN_CPU_CLUSTER_I_SIZE
#endif

/* Copies PBC shifted i-cell packed atom coordinates to working array */
static gmx_inline void
icell_set_x_simd_4xn(int ci,
                     real shx, real shy, real shz,
                     int gmx_unused stride, const real *x,
                     nbnxn_list_work_t *work)
{
    int    ia;
    real  *x_ci_simd = work->x_ci_simd;

    ia = xIndexFromCi<NbnxnLayout::Simd4xN>(ci);

    store(x_ci_simd +  0*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 0*STRIDE_S    ] + shx) );
    store(x_ci_simd +  1*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 1*STRIDE_S    ] + shy) );
    store(x_ci_simd +  2*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 2*STRIDE_S    ] + shz) );
    store(x_ci_simd +  3*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 0*STRIDE_S + 1] + shx) );
    store(x_ci_simd +  4*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 1*STRIDE_S + 1] + shy) );
    store(x_ci_simd +  5*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 2*STRIDE_S + 1] + shz) );
    store(x_ci_simd +  6*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 0*STRIDE_S + 2] + shx) );
    store(x_ci_simd +  7*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 1*STRIDE_S + 2] + shy) );
    store(x_ci_simd +  8*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 2*STRIDE_S + 2] + shz) );
    store(x_ci_simd +  9*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 0*STRIDE_S + 3] + shx) );
    store(x_ci_simd + 10*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 1*STRIDE_S + 3] + shy) );
    store(x_ci_simd + 11*GMX_SIMD_REAL_WIDTH, SimdReal(x[ia + 2*STRIDE_S + 3] + shz) );
}

/* SIMD code for checking and adding cluster-pairs to the list using coordinates in packed format.
 *
 * Checks bouding box distances and possibly atom pair distances.
 * This is an accelerated version of make_cluster_list_simple.
 *
 * \param[in]     gridj               The j-grid
 * \param[in,out] nbl                 The pair-list to store the cluster pairs in
 * \param[in]     icluster            The index of the i-cluster
 * \param[in]     firstCell           The first cluster in the j-range, using i-cluster size indexing
 * \param[in]     lastCell            The last cluster in the j-range, using i-cluster size indexing
 * \param[in]     excludeSubDiagonal  Exclude atom pairs with i-index > j-index
 * \param[in]     x_j                 Coordinates for the j-atom, in SIMD packed format
 * \param[in]     rlist2              The squared list cut-off
 * \param[in]     rbb2                The squared cut-off for putting cluster-pairs in the list based on bounding box distance only
 * \param[in,out] numDistanceChecks   The number of distance checks performed
 */
static gmx_inline void
makeClusterListSimd4xn(const nbnxn_grid_t *      gridj,
                       nbnxn_pairlist_t *        nbl,
                       int                       icluster,
                       int                       firstCell,
                       int                       lastCell,
                       bool                      excludeSubDiagonal,
                       const real * gmx_restrict x_j,
                       real                      rlist2,
                       float                     rbb2,
                       int * gmx_restrict        numDistanceChecks)
{
    const real * gmx_restrict          x_ci_simd = nbl->work->x_ci_simd;
    const nbnxn_bb_t * gmx_restrict    bb_ci     = nbl->work->bb_ci;

    SimdReal                           jx_S, jy_S, jz_S;

    SimdReal                           dx_S0, dy_S0, dz_S0;
    SimdReal                           dx_S1, dy_S1, dz_S1;
    SimdReal                           dx_S2, dy_S2, dz_S2;
    SimdReal                           dx_S3, dy_S3, dz_S3;

    SimdReal                           rsq_S0;
    SimdReal                           rsq_S1;
    SimdReal                           rsq_S2;
    SimdReal                           rsq_S3;

    SimdBool                           wco_S0;
    SimdBool                           wco_S1;
    SimdBool                           wco_S2;
    SimdBool                           wco_S3;
    SimdBool                           wco_any_S01, wco_any_S23, wco_any_S;

    SimdReal                           rc2_S;

    gmx_bool                           InRange;
    float                              d2;
    int                                xind_f, xind_l;

    /* Convert the j-range from i-cluster size indexing to j-cluster indexing */
    /* cppcheck-suppress selfAssignment . selfAssignment for width 4.*/
    int jclusterFirst = cjFromCi<NbnxnLayout::Simd4xN>(firstCell);
#if GMX_SIMD_REAL_WIDTH >= NBNXN_CPU_CLUSTER_I_SIZE
    int jclusterLast  = cjFromCi<NbnxnLayout::Simd4xN>(lastCell);
#else
    /* Set the correct last j-cluster with a j-cluster size of 2 */
    int jclusterLast  = cjFromCi<NbnxnLayout::Simd4xN>(lastCell + 1) - 1;
#endif
    GMX_ASSERT(jclusterLast >= jclusterFirst, "We should have a non-empty j-cluster range, since the calling code should have ensured a non-empty cell range");

    rc2_S   = SimdReal(rlist2);

    InRange = FALSE;
    while (!InRange && jclusterFirst <= jclusterLast)
    {
#ifdef NBNXN_SEARCH_BB_SIMD4
        d2 = subc_bb_dist2_simd4(0, bb_ci, jclusterFirst, gridj->bbj);
#else
        d2 = subc_bb_dist2(0, bb_ci, jclusterFirst, gridj->bbj);
#endif
        *numDistanceChecks += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rlist2)
        {
            xind_f  = xIndexFromCj<NbnxnLayout::Simd4xN>(cjFromCi<NbnxnLayout::Simd4xN>(gridj->cell0) + jclusterFirst);

            jx_S  = load(x_j + xind_f + 0*STRIDE_S);
            jy_S  = load(x_j + xind_f + 1*STRIDE_S);
            jz_S  = load(x_j + xind_f + 2*STRIDE_S);


            /* Calculate distance */
            dx_S0            = load(x_ci_simd +  0*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S0            = load(x_ci_simd +  1*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S0            = load(x_ci_simd +  2*GMX_SIMD_REAL_WIDTH) - jz_S;
            dx_S1            = load(x_ci_simd +  3*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S1            = load(x_ci_simd +  4*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S1            = load(x_ci_simd +  5*GMX_SIMD_REAL_WIDTH) - jz_S;
            dx_S2            = load(x_ci_simd +  6*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S2            = load(x_ci_simd +  7*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S2            = load(x_ci_simd +  8*GMX_SIMD_REAL_WIDTH) - jz_S;
            dx_S3            = load(x_ci_simd +  9*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S3            = load(x_ci_simd + 10*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S3            = load(x_ci_simd + 11*GMX_SIMD_REAL_WIDTH) - jz_S;

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_S0           = norm2(dx_S0, dy_S0, dz_S0);
            rsq_S1           = norm2(dx_S1, dy_S1, dz_S1);
            rsq_S2           = norm2(dx_S2, dy_S2, dz_S2);
            rsq_S3           = norm2(dx_S3, dy_S3, dz_S3);

            wco_S0           = (rsq_S0 < rc2_S);
            wco_S1           = (rsq_S1 < rc2_S);
            wco_S2           = (rsq_S2 < rc2_S);
            wco_S3           = (rsq_S3 < rc2_S);

            wco_any_S01      = wco_S0 || wco_S1;
            wco_any_S23      = wco_S2 || wco_S3;
            wco_any_S        = wco_any_S01 || wco_any_S23;

            InRange          = anyTrue(wco_any_S);

            *numDistanceChecks += 4*GMX_SIMD_REAL_WIDTH;
        }
        if (!InRange)
        {
            jclusterFirst++;
        }
    }
    if (!InRange)
    {
        return;
    }

    InRange = FALSE;
    while (!InRange && jclusterLast > jclusterFirst)
    {
#ifdef NBNXN_SEARCH_BB_SIMD4
        d2 = subc_bb_dist2_simd4(0, bb_ci, jclusterLast, gridj->bbj);
#else
        d2 = subc_bb_dist2(0, bb_ci, jclusterLast, gridj->bbj);
#endif
        *numDistanceChecks += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rlist2)
        {
            xind_l  = xIndexFromCj<NbnxnLayout::Simd4xN>(cjFromCi<NbnxnLayout::Simd4xN>(gridj->cell0) + jclusterLast);

            jx_S  = load(x_j +xind_l + 0*STRIDE_S);
            jy_S  = load(x_j +xind_l + 1*STRIDE_S);
            jz_S  = load(x_j +xind_l + 2*STRIDE_S);

            /* Calculate distance */
            dx_S0            = load(x_ci_simd +  0*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S0            = load(x_ci_simd +  1*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S0            = load(x_ci_simd +  2*GMX_SIMD_REAL_WIDTH) - jz_S;
            dx_S1            = load(x_ci_simd +  3*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S1            = load(x_ci_simd +  4*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S1            = load(x_ci_simd +  5*GMX_SIMD_REAL_WIDTH) - jz_S;
            dx_S2            = load(x_ci_simd +  6*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S2            = load(x_ci_simd +  7*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S2            = load(x_ci_simd +  8*GMX_SIMD_REAL_WIDTH) - jz_S;
            dx_S3            = load(x_ci_simd +  9*GMX_SIMD_REAL_WIDTH) - jx_S;
            dy_S3            = load(x_ci_simd + 10*GMX_SIMD_REAL_WIDTH) - jy_S;
            dz_S3            = load(x_ci_simd + 11*GMX_SIMD_REAL_WIDTH) - jz_S;

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_S0           = norm2(dx_S0, dy_S0, dz_S0);
            rsq_S1           = norm2(dx_S1, dy_S1, dz_S1);
            rsq_S2           = norm2(dx_S2, dy_S2, dz_S2);
            rsq_S3           = norm2(dx_S3, dy_S3, dz_S3);

            wco_S0           = (rsq_S0 < rc2_S);
            wco_S1           = (rsq_S1 < rc2_S);
            wco_S2           = (rsq_S2 < rc2_S);
            wco_S3           = (rsq_S3 < rc2_S);

            wco_any_S01      = wco_S0 || wco_S1;
            wco_any_S23      = wco_S2 || wco_S3;
            wco_any_S        = wco_any_S01 || wco_any_S23;

            InRange          = anyTrue(wco_any_S);

            *numDistanceChecks += 4*GMX_SIMD_REAL_WIDTH;
        }
        if (!InRange)
        {
            jclusterLast--;
        }
    }

    if (jclusterFirst <= jclusterLast)
    {
        for (int jcluster = jclusterFirst; jcluster <= jclusterLast; jcluster++)
        {
            /* Store cj and the interaction mask */
            nbl->cj[nbl->ncj].cj   = cjFromCi<NbnxnLayout::Simd4xN>(gridj->cell0) + jcluster;
            nbl->cj[nbl->ncj].excl = get_imask_simd_4xn(excludeSubDiagonal, icluster, jcluster);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

#undef STRIDE_S
