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

#include "gmxpre.h"

#include "nbnxn_grid.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_util.h"
#include "gromacs/utility/smalloc.h"

/* nbnxn_internal.h included gromacs/simd/macros.h */
#include "gromacs/mdlib/nbnxn_internal.h"
#ifdef GMX_SIMD
#include "gromacs/simd/vector_operations.h"
#endif

#ifdef NBNXN_SEARCH_BB_SIMD4
/* Always use 4-wide SIMD for bounding box calculations */

#    ifndef GMX_DOUBLE
/* Single precision BBs + coordinates, we can also load coordinates with SIMD */
#        define NBNXN_SEARCH_SIMD4_FLOAT_X_BB
#    endif

#    if defined NBNXN_SEARCH_SIMD4_FLOAT_X_BB && (GPU_NSUBCELL == 4 || GPU_NSUBCELL == 8)
/* Store bounding boxes with x, y and z coordinates in packs of 4 */
#        define NBNXN_PBB_SIMD4
#    endif

/* The packed bounding box coordinate stride is always set to 4.
 * With AVX we could use 8, but that turns out not to be faster.
 */
#    define STRIDE_PBB        4
#    define STRIDE_PBB_2LOG   2

#endif /* NBNXN_SEARCH_BB_SIMD4 */


#ifdef NBNXN_SEARCH_BB_SIMD4
/* Store bounding boxes corners as quadruplets: xxxxyyyyzzzz */
#define NBNXN_BBXXXX
/* Size of bounding box corners quadruplet */
#define NNBSBB_XXXX      (NNBSBB_D*DIM*STRIDE_PBB)
#endif


/* This define is a lazy way to avoid interdependence of the grid
 * and searching data structures.
 */
#define NBNXN_NA_SC_MAX (GPU_NSUBCELL*NBNXN_GPU_CLUSTER_SIZE)


static void nbnxn_grid_init(nbnxn_grid_t * grid)
{
    grid->cxy_na      = NULL;
    grid->cxy_ind     = NULL;
    grid->cxy_nalloc  = 0;
    grid->bb          = NULL;
    grid->bbj         = NULL;
    grid->nc_nalloc   = 0;
}

void nbnxn_grids_init(nbnxn_search_t nbs, int ngrid)
{
    int g;

    nbs->ngrid = ngrid;

    snew(nbs->grid, nbs->ngrid);
    for (g = 0; g < nbs->ngrid; g++)
    {
        nbnxn_grid_init(&nbs->grid[g]);
    }
}

static real grid_atom_density(int n, rvec corner0, rvec corner1)
{
    rvec size;

    if (n == 0)
    {
        /* To avoid zero density we use a minimum of 1 atom */
        n = 1;
    }

    rvec_sub(corner1, corner0, size);

    return n/(size[XX]*size[YY]*size[ZZ]);
}

/* Set the grid x/y cell dimensions and count for the local grid */
static int set_grid_size_xy(const nbnxn_search_t nbs,
                            nbnxn_grid_t *grid,
                            int n, rvec corner0, rvec corner1,
                            real atom_density)
{
    rvec size;
    int  na_c;
    real adens, tlen, tlen_x, tlen_y, nc_max;
    int  t;

    rvec_sub(corner1, corner0, size);

    if (n > grid->na_sc)
    {
        assert(atom_density > 0);

        /* target cell length */
        if (grid->bSimple)
        {
            /* To minimize the zero interactions, we should make
             * the largest of the i/j cell cubic.
             */
            na_c = max(grid->na_c, grid->na_cj);

            /* Approximately cubic cells */
            tlen   = pow(na_c/atom_density, 1.0/3.0);
            tlen_x = tlen;
            tlen_y = tlen;
        }
        else
        {
            /* Approximately cubic sub cells */
            tlen   = pow(grid->na_c/atom_density, 1.0/3.0);
            tlen_x = tlen*GPU_NSUBCELL_X;
            tlen_y = tlen*GPU_NSUBCELL_Y;
        }
        /* We don't round ncx and ncy to the nearest integer, but have
         * a preference for rouding down, because the pair search is cheaper
         * when the grid is coarser. Also there are somewhat fewer cell
         * pairs and filler particles in the list when the fixed cell
         * dimensions (x,y) are larger than the variable one (z), instead
         * of the other way around.
         * Note that the halo communication volume in x and y is on average
         * higher with coarser grids, but this depends strongly on the
         * exact values of the cut-off and the grid sizes.
         */
        grid->ncx = max(1, (int)(size[XX]/tlen_x + 0.25));
        grid->ncy = max(1, (int)(size[YY]/tlen_y + 0.25));
    }
    else
    {
        grid->ncx = 1;
        grid->ncy = 1;
    }

    grid->sx     = size[XX]/grid->ncx;
    grid->sy     = size[YY]/grid->ncy;
    grid->inv_sx = 1/grid->sx;
    grid->inv_sy = 1/grid->sy;

    if (debug)
    {
        fprintf(debug, "Local grid %d x %d cells of %.3f x %.3f nm\n",
                grid->ncx, grid->ncy, grid->sx, grid->sy);
    }

    /* We need one additional cell entry for particles moved by DD */
    if (grid->ncx*grid->ncy+1 > grid->cxy_nalloc)
    {
        grid->cxy_nalloc = over_alloc_large(grid->ncx*grid->ncy+1);
        srenew(grid->cxy_na, grid->cxy_nalloc);
        srenew(grid->cxy_ind, grid->cxy_nalloc+1);
    }
    for (t = 0; t < nbs->nthread_max; t++)
    {
        if (grid->ncx*grid->ncy+1 > nbs->work[t].cxy_na_nalloc)
        {
            nbs->work[t].cxy_na_nalloc = over_alloc_large(grid->ncx*grid->ncy+1);
            srenew(nbs->work[t].cxy_na, nbs->work[t].cxy_na_nalloc);
        }
    }

    /* Worst case scenario of 1 atom in each last cell */
    if (grid->na_cj <= grid->na_c)
    {
        nc_max = n/grid->na_sc + grid->ncx*grid->ncy;
    }
    else
    {
        nc_max = n/grid->na_sc + grid->ncx*grid->ncy*grid->na_cj/grid->na_c;
    }

    copy_rvec(corner0, grid->c0);
    copy_rvec(corner1, grid->c1);
    copy_rvec(size,    grid->size);

    return nc_max;
}

/* Reallocate bounding boxes and flags for grid, when necessary */
static void grid_realloc_bb_flags(nbnxn_grid_t *grid, gmx_bool bFEP, int nc_max)
{
    if (nc_max > grid->nc_nalloc)
    {
        grid->nc_nalloc = over_alloc_large(nc_max);
        srenew(grid->nsubc, grid->nc_nalloc);
        srenew(grid->bbcz, grid->nc_nalloc*NNBSBB_D);

        sfree_aligned(grid->bb);
        /* This snew also zeros the contents, this avoid possible
         * floating exceptions in SIMD with the unused bb elements.
         */
        if (grid->bSimple)
        {
            snew_aligned(grid->bb, grid->nc_nalloc, 16);
        }
        else
        {
#ifdef NBNXN_BBXXXX
            int pbb_nalloc;

            pbb_nalloc = grid->nc_nalloc*GPU_NSUBCELL/STRIDE_PBB*NNBSBB_XXXX;
            snew_aligned(grid->pbb, pbb_nalloc, 16);
#else
            snew_aligned(grid->bb, grid->nc_nalloc*GPU_NSUBCELL, 16);
#endif
        }

        if (grid->bSimple)
        {
            if (grid->na_cj == grid->na_c)
            {
                grid->bbj = grid->bb;
            }
            else
            {
                sfree_aligned(grid->bbj);
                snew_aligned(grid->bbj, grid->nc_nalloc*grid->na_c/grid->na_cj, 16);
            }
        }

        srenew(grid->flags, grid->nc_nalloc);
        if (bFEP)
        {
            srenew(grid->fep, grid->nc_nalloc*grid->na_sc/grid->na_c);
        }
    }
}

/* We need to sort paricles in grid columns on z-coordinate.
 * As particle are very often distributed homogeneously, we a sorting
 * algorithm similar to pigeonhole sort. We multiply the z-coordinate
 * by a factor, cast to an int and try to store in that hole. If the hole
 * is full, we move this or another particle. A second pass is needed to make
 * contiguous elements. SORT_GRID_OVERSIZE is the ratio of holes to particles.
 * 4 is the optimal value for homogeneous particle distribution and allows
 * for an O(#particles) sort up till distributions were all particles are
 * concentrated in 1/4 of the space. No NlogN fallback is implemented,
 * as it can be expensive to detect imhomogeneous particle distributions.
 * SGSF is the maximum ratio of holes used, in the worst case all particles
 * end up in the last hole and we need #particles extra holes at the end.
 */
#define SORT_GRID_OVERSIZE 4
#define SGSF (SORT_GRID_OVERSIZE + 1)

/* Sort particle index a on coordinates x along dim.
 * Backwards tells if we want decreasing iso increasing coordinates.
 * h0 is the minimum of the coordinate range.
 * invh is the 1/length of the sorting range.
 * n_per_h (>=n) is the expected average number of particles per 1/invh
 * sort is the sorting work array.
 * sort should have a size of at least n_per_h*SORT_GRID_OVERSIZE + n,
 * or easier, allocate at least n*SGSF elements.
 */
static void sort_atoms(int dim, gmx_bool Backwards,
                       int gmx_unused dd_zone,
                       int *a, int n, rvec *x,
                       real h0, real invh, int n_per_h,
                       int *sort)
{
    int nsort, i, c;
    int zi, zim, zi_min, zi_max;
    int cp, tmp;

    if (n <= 1)
    {
        /* Nothing to do */
        return;
    }

#ifndef NDEBUG
    if (n > n_per_h)
    {
        gmx_incons("n > n_per_h");
    }
#endif

    /* Transform the inverse range height into the inverse hole height */
    invh *= n_per_h*SORT_GRID_OVERSIZE;

    /* Set nsort to the maximum possible number of holes used.
     * In worst case all n elements end up in the last bin.
     */
    nsort = n_per_h*SORT_GRID_OVERSIZE + n;

    /* Determine the index range used, so we can limit it for the second pass */
    zi_min = INT_MAX;
    zi_max = -1;

    /* Sort the particles using a simple index sort */
    for (i = 0; i < n; i++)
    {
        /* The cast takes care of float-point rounding effects below zero.
         * This code assumes particles are less than 1/SORT_GRID_OVERSIZE
         * times the box height out of the box.
         */
        zi = (int)((x[a[i]][dim] - h0)*invh);

#ifndef NDEBUG
        /* As we can have rounding effect, we use > iso >= here */
        if (zi < 0 || (dd_zone == 0 && zi > n_per_h*SORT_GRID_OVERSIZE))
        {
            gmx_fatal(FARGS, "(int)((x[%d][%c]=%f - %f)*%f) = %d, not in 0 - %d*%d\n",
                      a[i], 'x'+dim, x[a[i]][dim], h0, invh, zi,
                      n_per_h, SORT_GRID_OVERSIZE);
        }
#endif

        /* In a non-local domain, particles communcated for bonded interactions
         * can be far beyond the grid size, which is set by the non-bonded
         * cut-off distance. We sort such particles into the last cell.
         */
        if (zi > n_per_h*SORT_GRID_OVERSIZE)
        {
            zi = n_per_h*SORT_GRID_OVERSIZE;
        }

        /* Ideally this particle should go in sort cell zi,
         * but that might already be in use,
         * in that case find the first empty cell higher up
         */
        if (sort[zi] < 0)
        {
            sort[zi] = a[i];
            zi_min   = min(zi_min, zi);
            zi_max   = max(zi_max, zi);
        }
        else
        {
            /* We have multiple atoms in the same sorting slot.
             * Sort on real z for minimal bounding box size.
             * There is an extra check for identical z to ensure
             * well-defined output order, independent of input order
             * to ensure binary reproducibility after restarts.
             */
            while (sort[zi] >= 0 && ( x[a[i]][dim] >  x[sort[zi]][dim] ||
                                      (x[a[i]][dim] == x[sort[zi]][dim] &&
                                       a[i] > sort[zi])))
            {
                zi++;
            }

            if (sort[zi] >= 0)
            {
                /* Shift all elements by one slot until we find an empty slot */
                cp  = sort[zi];
                zim = zi + 1;
                while (sort[zim] >= 0)
                {
                    tmp       = sort[zim];
                    sort[zim] = cp;
                    cp        = tmp;
                    zim++;
                }
                sort[zim] = cp;
                zi_max    = max(zi_max, zim);
            }
            sort[zi] = a[i];
            zi_max   = max(zi_max, zi);
        }
    }

    c = 0;
    if (!Backwards)
    {
        for (zi = 0; zi < nsort; zi++)
        {
            if (sort[zi] >= 0)
            {
                a[c++]   = sort[zi];
                sort[zi] = -1;
            }
        }
    }
    else
    {
        for (zi = zi_max; zi >= zi_min; zi--)
        {
            if (sort[zi] >= 0)
            {
                a[c++]   = sort[zi];
                sort[zi] = -1;
            }
        }
    }
    if (c < n)
    {
        gmx_incons("Lost particles while sorting");
    }
}

#ifdef GMX_DOUBLE
#define R2F_D(x) ((float)((x) >= 0 ? ((1-GMX_FLOAT_EPS)*(x)) : ((1+GMX_FLOAT_EPS)*(x))))
#define R2F_U(x) ((float)((x) >= 0 ? ((1+GMX_FLOAT_EPS)*(x)) : ((1-GMX_FLOAT_EPS)*(x))))
#else
#define R2F_D(x) (x)
#define R2F_U(x) (x)
#endif

/* Coordinate order x,y,z, bb order xyz0 */
static void calc_bounding_box(int na, int stride, const real *x, nbnxn_bb_t *bb)
{
    int  i, j;
    real xl, xh, yl, yh, zl, zh;

    i  = 0;
    xl = x[i+XX];
    xh = x[i+XX];
    yl = x[i+YY];
    yh = x[i+YY];
    zl = x[i+ZZ];
    zh = x[i+ZZ];
    i += stride;
    for (j = 1; j < na; j++)
    {
        xl = min(xl, x[i+XX]);
        xh = max(xh, x[i+XX]);
        yl = min(yl, x[i+YY]);
        yh = max(yh, x[i+YY]);
        zl = min(zl, x[i+ZZ]);
        zh = max(zh, x[i+ZZ]);
        i += stride;
    }
    /* Note: possible double to float conversion here */
    bb->lower[BB_X] = R2F_D(xl);
    bb->lower[BB_Y] = R2F_D(yl);
    bb->lower[BB_Z] = R2F_D(zl);
    bb->upper[BB_X] = R2F_U(xh);
    bb->upper[BB_Y] = R2F_U(yh);
    bb->upper[BB_Z] = R2F_U(zh);
}

/* Packed coordinates, bb order xyz0 */
static void calc_bounding_box_x_x4(int na, const real *x, nbnxn_bb_t *bb)
{
    int  j;
    real xl, xh, yl, yh, zl, zh;

    xl = x[XX*PACK_X4];
    xh = x[XX*PACK_X4];
    yl = x[YY*PACK_X4];
    yh = x[YY*PACK_X4];
    zl = x[ZZ*PACK_X4];
    zh = x[ZZ*PACK_X4];
    for (j = 1; j < na; j++)
    {
        xl = min(xl, x[j+XX*PACK_X4]);
        xh = max(xh, x[j+XX*PACK_X4]);
        yl = min(yl, x[j+YY*PACK_X4]);
        yh = max(yh, x[j+YY*PACK_X4]);
        zl = min(zl, x[j+ZZ*PACK_X4]);
        zh = max(zh, x[j+ZZ*PACK_X4]);
    }
    /* Note: possible double to float conversion here */
    bb->lower[BB_X] = R2F_D(xl);
    bb->lower[BB_Y] = R2F_D(yl);
    bb->lower[BB_Z] = R2F_D(zl);
    bb->upper[BB_X] = R2F_U(xh);
    bb->upper[BB_Y] = R2F_U(yh);
    bb->upper[BB_Z] = R2F_U(zh);
}

/* Packed coordinates, bb order xyz0 */
static void calc_bounding_box_x_x8(int na, const real *x, nbnxn_bb_t *bb)
{
    int  j;
    real xl, xh, yl, yh, zl, zh;

    xl = x[XX*PACK_X8];
    xh = x[XX*PACK_X8];
    yl = x[YY*PACK_X8];
    yh = x[YY*PACK_X8];
    zl = x[ZZ*PACK_X8];
    zh = x[ZZ*PACK_X8];
    for (j = 1; j < na; j++)
    {
        xl = min(xl, x[j+XX*PACK_X8]);
        xh = max(xh, x[j+XX*PACK_X8]);
        yl = min(yl, x[j+YY*PACK_X8]);
        yh = max(yh, x[j+YY*PACK_X8]);
        zl = min(zl, x[j+ZZ*PACK_X8]);
        zh = max(zh, x[j+ZZ*PACK_X8]);
    }
    /* Note: possible double to float conversion here */
    bb->lower[BB_X] = R2F_D(xl);
    bb->lower[BB_Y] = R2F_D(yl);
    bb->lower[BB_Z] = R2F_D(zl);
    bb->upper[BB_X] = R2F_U(xh);
    bb->upper[BB_Y] = R2F_U(yh);
    bb->upper[BB_Z] = R2F_U(zh);
}

/* Packed coordinates, bb order xyz0 */
static void calc_bounding_box_x_x4_halves(int na, const real *x,
                                          nbnxn_bb_t *bb, nbnxn_bb_t *bbj)
{
    calc_bounding_box_x_x4(min(na, 2), x, bbj);

    if (na > 2)
    {
        calc_bounding_box_x_x4(min(na-2, 2), x+(PACK_X4>>1), bbj+1);
    }
    else
    {
        /* Set the "empty" bounding box to the same as the first one,
         * so we don't need to treat special cases in the rest of the code.
         */
#ifdef NBNXN_SEARCH_BB_SIMD4
        gmx_simd4_store_f(&bbj[1].lower[0], gmx_simd4_load_f(&bbj[0].lower[0]));
        gmx_simd4_store_f(&bbj[1].upper[0], gmx_simd4_load_f(&bbj[0].upper[0]));
#else
        bbj[1] = bbj[0];
#endif
    }

#ifdef NBNXN_SEARCH_BB_SIMD4
    gmx_simd4_store_f(&bb->lower[0],
                      gmx_simd4_min_f(gmx_simd4_load_f(&bbj[0].lower[0]),
                                      gmx_simd4_load_f(&bbj[1].lower[0])));
    gmx_simd4_store_f(&bb->upper[0],
                      gmx_simd4_max_f(gmx_simd4_load_f(&bbj[0].upper[0]),
                                      gmx_simd4_load_f(&bbj[1].upper[0])));
#else
    {
        int i;

        for (i = 0; i < NNBSBB_C; i++)
        {
            bb->lower[i] = min(bbj[0].lower[i], bbj[1].lower[i]);
            bb->upper[i] = max(bbj[0].upper[i], bbj[1].upper[i]);
        }
    }
#endif
}

#ifdef NBNXN_SEARCH_BB_SIMD4

/* Coordinate order xyz, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx(int na, int stride, const real *x, float *bb)
{
    int  i, j;
    real xl, xh, yl, yh, zl, zh;

    i  = 0;
    xl = x[i+XX];
    xh = x[i+XX];
    yl = x[i+YY];
    yh = x[i+YY];
    zl = x[i+ZZ];
    zh = x[i+ZZ];
    i += stride;
    for (j = 1; j < na; j++)
    {
        xl = min(xl, x[i+XX]);
        xh = max(xh, x[i+XX]);
        yl = min(yl, x[i+YY]);
        yh = max(yh, x[i+YY]);
        zl = min(zl, x[i+ZZ]);
        zh = max(zh, x[i+ZZ]);
        i += stride;
    }
    /* Note: possible double to float conversion here */
    bb[0*STRIDE_PBB] = R2F_D(xl);
    bb[1*STRIDE_PBB] = R2F_D(yl);
    bb[2*STRIDE_PBB] = R2F_D(zl);
    bb[3*STRIDE_PBB] = R2F_U(xh);
    bb[4*STRIDE_PBB] = R2F_U(yh);
    bb[5*STRIDE_PBB] = R2F_U(zh);
}

#endif /* NBNXN_SEARCH_BB_SIMD4 */

#ifdef NBNXN_SEARCH_SIMD4_FLOAT_X_BB

/* Coordinate order xyz?, bb order xyz0 */
static void calc_bounding_box_simd4(int na, const float *x, nbnxn_bb_t *bb)
{
    gmx_simd4_float_t bb_0_S, bb_1_S;
    gmx_simd4_float_t x_S;

    int               i;

    bb_0_S = gmx_simd4_load_f(x);
    bb_1_S = bb_0_S;

    for (i = 1; i < na; i++)
    {
        x_S    = gmx_simd4_load_f(x+i*NNBSBB_C);
        bb_0_S = gmx_simd4_min_f(bb_0_S, x_S);
        bb_1_S = gmx_simd4_max_f(bb_1_S, x_S);
    }

    gmx_simd4_store_f(&bb->lower[0], bb_0_S);
    gmx_simd4_store_f(&bb->upper[0], bb_1_S);
}

/* Coordinate order xyz?, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx_simd4(int na, const float *x,
                                         nbnxn_bb_t *bb_work_aligned,
                                         real *bb)
{
    calc_bounding_box_simd4(na, x, bb_work_aligned);

    bb[0*STRIDE_PBB] = bb_work_aligned->lower[BB_X];
    bb[1*STRIDE_PBB] = bb_work_aligned->lower[BB_Y];
    bb[2*STRIDE_PBB] = bb_work_aligned->lower[BB_Z];
    bb[3*STRIDE_PBB] = bb_work_aligned->upper[BB_X];
    bb[4*STRIDE_PBB] = bb_work_aligned->upper[BB_Y];
    bb[5*STRIDE_PBB] = bb_work_aligned->upper[BB_Z];
}

#endif /* NBNXN_SEARCH_SIMD4_FLOAT_X_BB */


/* Combines pairs of consecutive bounding boxes */
static void combine_bounding_box_pairs(nbnxn_grid_t *grid, const nbnxn_bb_t *bb)
{
    int    i, j, sc2, nc2, c2;

    for (i = 0; i < grid->ncx*grid->ncy; i++)
    {
        /* Starting bb in a column is expected to be 2-aligned */
        sc2 = grid->cxy_ind[i]>>1;
        /* For odd numbers skip the last bb here */
        nc2 = (grid->cxy_na[i]+3)>>(2+1);
        for (c2 = sc2; c2 < sc2+nc2; c2++)
        {
#ifdef NBNXN_SEARCH_BB_SIMD4
            gmx_simd4_float_t min_S, max_S;

            min_S = gmx_simd4_min_f(gmx_simd4_load_f(&bb[c2*2+0].lower[0]),
                                    gmx_simd4_load_f(&bb[c2*2+1].lower[0]));
            max_S = gmx_simd4_max_f(gmx_simd4_load_f(&bb[c2*2+0].upper[0]),
                                    gmx_simd4_load_f(&bb[c2*2+1].upper[0]));
            gmx_simd4_store_f(&grid->bbj[c2].lower[0], min_S);
            gmx_simd4_store_f(&grid->bbj[c2].upper[0], max_S);
#else
            for (j = 0; j < NNBSBB_C; j++)
            {
                grid->bbj[c2].lower[j] = min(bb[c2*2+0].lower[j],
                                             bb[c2*2+1].lower[j]);
                grid->bbj[c2].upper[j] = max(bb[c2*2+0].upper[j],
                                             bb[c2*2+1].upper[j]);
            }
#endif
        }
        if (((grid->cxy_na[i]+3)>>2) & 1)
        {
            /* The bb count in this column is odd: duplicate the last bb */
            for (j = 0; j < NNBSBB_C; j++)
            {
                grid->bbj[c2].lower[j] = bb[c2*2].lower[j];
                grid->bbj[c2].upper[j] = bb[c2*2].upper[j];
            }
        }
    }
}


/* Prints the average bb size, used for debug output */
static void print_bbsizes_simple(FILE                *fp,
                                 const nbnxn_grid_t  *grid)
{
    int  c, d;
    dvec ba;

    clear_dvec(ba);
    for (c = 0; c < grid->nc; c++)
    {
        for (d = 0; d < DIM; d++)
        {
            ba[d] += grid->bb[c].upper[d] - grid->bb[c].lower[d];
        }
    }
    dsvmul(1.0/grid->nc, ba, ba);

    fprintf(fp, "ns bb: grid %4.2f %4.2f %4.2f abs %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            grid->sx,
            grid->sy,
            grid->na_c/(grid->atom_density*grid->sx*grid->sy),
            ba[XX], ba[YY], ba[ZZ],
            ba[XX]/grid->sx,
            ba[YY]/grid->sy,
            ba[ZZ]/(grid->na_c/(grid->atom_density*grid->sx*grid->sy)));
}

/* Prints the average bb size, used for debug output */
static void print_bbsizes_supersub(FILE                *fp,
                                   const nbnxn_grid_t  *grid)
{
    int  ns, c, s;
    dvec ba;

    clear_dvec(ba);
    ns = 0;
    for (c = 0; c < grid->nc; c++)
    {
#ifdef NBNXN_BBXXXX
        for (s = 0; s < grid->nsubc[c]; s += STRIDE_PBB)
        {
            int cs_w, i, d;

            cs_w = (c*GPU_NSUBCELL + s)/STRIDE_PBB;
            for (i = 0; i < STRIDE_PBB; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    ba[d] +=
                        grid->pbb[cs_w*NNBSBB_XXXX+(DIM+d)*STRIDE_PBB+i] -
                        grid->pbb[cs_w*NNBSBB_XXXX+     d *STRIDE_PBB+i];
                }
            }
        }
#else
        for (s = 0; s < grid->nsubc[c]; s++)
        {
            int cs, d;

            cs = c*GPU_NSUBCELL + s;
            for (d = 0; d < DIM; d++)
            {
                ba[d] += grid->bb[cs].upper[d] - grid->bb[cs].lower[d];
            }
        }
#endif
        ns += grid->nsubc[c];
    }
    dsvmul(1.0/ns, ba, ba);

    fprintf(fp, "ns bb: grid %4.2f %4.2f %4.2f abs %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            grid->sx/GPU_NSUBCELL_X,
            grid->sy/GPU_NSUBCELL_Y,
            grid->na_sc/(grid->atom_density*grid->sx*grid->sy*GPU_NSUBCELL_Z),
            ba[XX], ba[YY], ba[ZZ],
            ba[XX]*GPU_NSUBCELL_X/grid->sx,
            ba[YY]*GPU_NSUBCELL_Y/grid->sy,
            ba[ZZ]/(grid->na_sc/(grid->atom_density*grid->sx*grid->sy*GPU_NSUBCELL_Z)));
}

/* Set non-bonded interaction flags for the current cluster.
 * Sorts atoms on LJ coefficients: !=0 first, ==0 at the end.
 */
static void sort_cluster_on_flag(int na_c,
                                 int a0, int a1, const int *atinfo,
                                 int *order,
                                 int *flags)
{
    int      subc, s, a, n1, n2, a_lj_max, i, j;
    int      sort1[NBNXN_NA_SC_MAX/GPU_NSUBCELL];
    int      sort2[NBNXN_NA_SC_MAX/GPU_NSUBCELL];
    gmx_bool haveQ, bFEP;

    *flags = 0;

    subc = 0;
    for (s = a0; s < a1; s += na_c)
    {
        /* Make lists for this (sub-)cell on atoms with and without LJ */
        n1       = 0;
        n2       = 0;
        haveQ    = FALSE;
        a_lj_max = -1;
        for (a = s; a < min(s+na_c, a1); a++)
        {
            haveQ = haveQ || GET_CGINFO_HAS_Q(atinfo[order[a]]);

            if (GET_CGINFO_HAS_VDW(atinfo[order[a]]))
            {
                sort1[n1++] = order[a];
                a_lj_max    = a;
            }
            else
            {
                sort2[n2++] = order[a];
            }
        }

        /* If we don't have atoms with LJ, there's nothing to sort */
        if (n1 > 0)
        {
            *flags |= NBNXN_CI_DO_LJ(subc);

            if (2*n1 <= na_c)
            {
                /* Only sort when strictly necessary.
                 * Ordering particles can lead to less accurate summation
                 * due to rounding, both for LJ and Coulomb interactions.
                 */
                if (2*(a_lj_max - s) >= na_c)
                {
                    for (i = 0; i < n1; i++)
                    {
                        order[a0+i] = sort1[i];
                    }
                    for (j = 0; j < n2; j++)
                    {
                        order[a0+n1+j] = sort2[j];
                    }
                }

                *flags |= NBNXN_CI_HALF_LJ(subc);
            }
        }
        if (haveQ)
        {
            *flags |= NBNXN_CI_DO_COUL(subc);
        }
        subc++;
    }
}

/* Fill a pair search cell with atoms.
 * Potentially sorts atoms and sets the interaction flags.
 */
void fill_cell(const nbnxn_search_t nbs,
               nbnxn_grid_t *grid,
               nbnxn_atomdata_t *nbat,
               int a0, int a1,
               const int *atinfo,
               rvec *x,
               int sx, int sy, int sz,
               nbnxn_bb_t gmx_unused *bb_work_aligned)
{
    int         na, a;
    size_t      offset;
    nbnxn_bb_t *bb_ptr;
#ifdef NBNXN_BBXXXX
    float      *pbb_ptr;
#endif

    na = a1 - a0;

    if (grid->bSimple)
    {
        /* Note that non-local grids are already sorted. Then flag_sort_cluster
         * will only set the flags and no actual sorting will happen.
         */
        sort_cluster_on_flag(grid->na_c, a0, a1, atinfo, nbs->a,
                             grid->flags+(a0>>grid->na_c_2log)-grid->cell0);
    }

    if (nbs->bFEP)
    {
        /* Set the fep flag for perturbed atoms in this (sub-)cell */
        int c, at;

        /* The grid-local cluster/(sub-)cell index */
        c            = (a0 >> grid->na_c_2log) - grid->cell0*(grid->bSimple ? 1 : GPU_NSUBCELL);
        grid->fep[c] = 0;
        for (at = a0; at < a1; at++)
        {
            if (nbs->a[at] >= 0 && GET_CGINFO_FEP(atinfo[nbs->a[at]]))
            {
                grid->fep[c] |= (1 << (at - a0));
            }
        }
    }

    /* Now we have sorted the atoms, set the cell indices */
    for (a = a0; a < a1; a++)
    {
        nbs->cell[nbs->a[a]] = a;
    }

    copy_rvec_to_nbat_real(nbs->a+a0, a1-a0, grid->na_c, x,
                           nbat->XFormat, nbat->x, a0,
                           sx, sy, sz);

    if (nbat->XFormat == nbatX4)
    {
        /* Store the bounding boxes as xyz.xyz. */
        offset = (a0 - grid->cell0*grid->na_sc) >> grid->na_c_2log;
        bb_ptr = grid->bb + offset;

#if defined GMX_NBNXN_SIMD && GMX_SIMD_REAL_WIDTH == 2
        if (2*grid->na_cj == grid->na_c)
        {
            calc_bounding_box_x_x4_halves(na, nbat->x+X4_IND_A(a0), bb_ptr,
                                          grid->bbj+offset*2);
        }
        else
#endif
        {
            calc_bounding_box_x_x4(na, nbat->x+X4_IND_A(a0), bb_ptr);
        }
    }
    else if (nbat->XFormat == nbatX8)
    {
        /* Store the bounding boxes as xyz.xyz. */
        offset = (a0 - grid->cell0*grid->na_sc) >> grid->na_c_2log;
        bb_ptr = grid->bb + offset;

        calc_bounding_box_x_x8(na, nbat->x+X8_IND_A(a0), bb_ptr);
    }
#ifdef NBNXN_BBXXXX
    else if (!grid->bSimple)
    {
        /* Store the bounding boxes in a format convenient
         * for SIMD4 calculations: xxxxyyyyzzzz...
         */
        pbb_ptr =
            grid->pbb +
            ((a0-grid->cell0*grid->na_sc)>>(grid->na_c_2log+STRIDE_PBB_2LOG))*NNBSBB_XXXX +
            (((a0-grid->cell0*grid->na_sc)>>grid->na_c_2log) & (STRIDE_PBB-1));

#ifdef NBNXN_SEARCH_SIMD4_FLOAT_X_BB
        if (nbat->XFormat == nbatXYZQ)
        {
            calc_bounding_box_xxxx_simd4(na, nbat->x+a0*nbat->xstride,
                                         bb_work_aligned, pbb_ptr);
        }
        else
#endif
        {
            calc_bounding_box_xxxx(na, nbat->xstride, nbat->x+a0*nbat->xstride,
                                   pbb_ptr);
        }
        if (gmx_debug_at)
        {
            fprintf(debug, "%2d %2d %2d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    sx, sy, sz,
                    pbb_ptr[0*STRIDE_PBB], pbb_ptr[3*STRIDE_PBB],
                    pbb_ptr[1*STRIDE_PBB], pbb_ptr[4*STRIDE_PBB],
                    pbb_ptr[2*STRIDE_PBB], pbb_ptr[5*STRIDE_PBB]);
        }
    }
#endif
    else
    {
        /* Store the bounding boxes as xyz.xyz. */
        bb_ptr = grid->bb+((a0-grid->cell0*grid->na_sc)>>grid->na_c_2log);

        calc_bounding_box(na, nbat->xstride, nbat->x+a0*nbat->xstride,
                          bb_ptr);

        if (gmx_debug_at)
        {
            int bbo;
            bbo = (a0 - grid->cell0*grid->na_sc)/grid->na_c;
            fprintf(debug, "%2d %2d %2d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    sx, sy, sz,
                    grid->bb[bbo].lower[BB_X],
                    grid->bb[bbo].lower[BB_Y],
                    grid->bb[bbo].lower[BB_Z],
                    grid->bb[bbo].upper[BB_X],
                    grid->bb[bbo].upper[BB_Y],
                    grid->bb[bbo].upper[BB_Z]);
        }
    }
}

/* Spatially sort the atoms within one grid column */
static void sort_columns_simple(const nbnxn_search_t nbs,
                                int dd_zone,
                                nbnxn_grid_t *grid,
                                int a0, int a1,
                                const int *atinfo,
                                rvec *x,
                                nbnxn_atomdata_t *nbat,
                                int cxy_start, int cxy_end,
                                int *sort_work)
{
    int  cxy;
    int  cx, cy, cz, ncz, cfilled, c;
    int  na, ash, ind, a;
    int  na_c, ash_c;

    if (debug)
    {
        fprintf(debug, "cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0, cxy_start, cxy_end, a0, a1);
    }

    /* Sort the atoms within each x,y column in 3 dimensions */
    for (cxy = cxy_start; cxy < cxy_end; cxy++)
    {
        cx = cxy/grid->ncy;
        cy = cxy - cx*grid->ncy;

        na  = grid->cxy_na[cxy];
        ncz = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

        if (dd_zone == 0)
        {
            /* Sort the atoms within each x,y column on z coordinate */
            sort_atoms(ZZ, FALSE, dd_zone,
                       nbs->a+ash, na, x,
                       grid->c0[ZZ],
                       1.0/grid->size[ZZ], ncz*grid->na_sc,
                       sort_work);
        }

        /* Fill the ncz cells in this column */
        cfilled = grid->cxy_ind[cxy];
        for (cz = 0; cz < ncz; cz++)
        {
            c  = grid->cxy_ind[cxy] + cz;

            ash_c = ash + cz*grid->na_sc;
            na_c  = min(grid->na_sc, na-(ash_c-ash));

            fill_cell(nbs, grid, nbat,
                      ash_c, ash_c+na_c, atinfo, x,
                      grid->na_sc*cx + (dd_zone >> 2),
                      grid->na_sc*cy + (dd_zone & 3),
                      grid->na_sc*cz,
                      NULL);

            /* This copy to bbcz is not really necessary.
             * But it allows to use the same grid search code
             * for the simple and supersub cell setups.
             */
            if (na_c > 0)
            {
                cfilled = c;
            }
            grid->bbcz[c*NNBSBB_D  ] = grid->bb[cfilled].lower[BB_Z];
            grid->bbcz[c*NNBSBB_D+1] = grid->bb[cfilled].upper[BB_Z];
        }

        /* Set the unused atom indices to -1 */
        for (ind = na; ind < ncz*grid->na_sc; ind++)
        {
            nbs->a[ash+ind] = -1;
        }
    }
}

/* Spatially sort the atoms within one grid column */
static void sort_columns_supersub(const nbnxn_search_t nbs,
                                  int dd_zone,
                                  nbnxn_grid_t *grid,
                                  int a0, int a1,
                                  const int *atinfo,
                                  rvec *x,
                                  nbnxn_atomdata_t *nbat,
                                  int cxy_start, int cxy_end,
                                  int *sort_work)
{
    int        cxy;
    int        cx, cy, cz = -1, c = -1, ncz;
    int        na, ash, na_c, ind, a;
    int        subdiv_z, sub_z, na_z, ash_z;
    int        subdiv_y, sub_y, na_y, ash_y;
    int        subdiv_x, sub_x, na_x, ash_x;

    nbnxn_bb_t bb_work_array[2], *bb_work_aligned;

    bb_work_aligned = (nbnxn_bb_t *)(((size_t)(bb_work_array+1)) & (~((size_t)15)));

    if (debug)
    {
        fprintf(debug, "cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0, cxy_start, cxy_end, a0, a1);
    }

    subdiv_x = grid->na_c;
    subdiv_y = GPU_NSUBCELL_X*subdiv_x;
    subdiv_z = GPU_NSUBCELL_Y*subdiv_y;

    /* Sort the atoms within each x,y column in 3 dimensions */
    for (cxy = cxy_start; cxy < cxy_end; cxy++)
    {
        cx = cxy/grid->ncy;
        cy = cxy - cx*grid->ncy;

        na  = grid->cxy_na[cxy];
        ncz = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

        if (dd_zone == 0)
        {
            /* Sort the atoms within each x,y column on z coordinate */
            sort_atoms(ZZ, FALSE, dd_zone,
                       nbs->a+ash, na, x,
                       grid->c0[ZZ],
                       1.0/grid->size[ZZ], ncz*grid->na_sc,
                       sort_work);
        }

        /* This loop goes over the supercells and subcells along z at once */
        for (sub_z = 0; sub_z < ncz*GPU_NSUBCELL_Z; sub_z++)
        {
            ash_z = ash + sub_z*subdiv_z;
            na_z  = min(subdiv_z, na-(ash_z-ash));

            /* We have already sorted on z */

            if (sub_z % GPU_NSUBCELL_Z == 0)
            {
                cz = sub_z/GPU_NSUBCELL_Z;
                c  = grid->cxy_ind[cxy] + cz;

                /* The number of atoms in this supercell */
                na_c = min(grid->na_sc, na-(ash_z-ash));

                grid->nsubc[c] = min(GPU_NSUBCELL, (na_c+grid->na_c-1)/grid->na_c);

                /* Store the z-boundaries of the super cell */
                if (dd_zone == 0)
                {
                    /* The particles are now sorted on z: use first and last */
                    grid->bbcz[c*NNBSBB_D  ] = x[nbs->a[ash_z]][ZZ];
                    grid->bbcz[c*NNBSBB_D+1] = x[nbs->a[ash_z+na_c-1]][ZZ];
                }
                else
                {
                    /* Within a cell the particles are not sorted on z */
                    real lower, upper;
                    int  i;

                    lower = x[nbs->a[ash_z]][ZZ];
                    upper = x[nbs->a[ash_z]][ZZ];
                    for (i = ash_z + 1; i < ash_z + na_c; i++)
                    {
                        lower = min(lower, x[nbs->a[i]][ZZ]);
                        upper = max(upper, x[nbs->a[i]][ZZ]);
                    }
                    grid->bbcz[c*NNBSBB_D  ] = lower;
                    grid->bbcz[c*NNBSBB_D+1] = upper;
                }
            }

#if GPU_NSUBCELL_Y > 1
            if (dd_zone == 0)
            {
                /* Sort the atoms along y */
                sort_atoms(YY, (sub_z & 1), dd_zone,
                           nbs->a+ash_z, na_z, x,
                           grid->c0[YY]+cy*grid->sy,
                           grid->inv_sy, subdiv_z,
                           sort_work);
            }
#endif

            for (sub_y = 0; sub_y < GPU_NSUBCELL_Y; sub_y++)
            {
                ash_y = ash_z + sub_y*subdiv_y;
                na_y  = min(subdiv_y, na-(ash_y-ash));

#if GPU_NSUBCELL_X > 1
                if (dd_zone == 0)
                {
                    /* Sort the atoms along x */
                    sort_atoms(XX, ((cz*GPU_NSUBCELL_Y + sub_y) & 1), dd_zone,
                               nbs->a+ash_y, na_y, x,
                               grid->c0[XX]+cx*grid->sx,
                               grid->inv_sx, subdiv_y,
                               sort_work);
                }
#endif

                for (sub_x = 0; sub_x < GPU_NSUBCELL_X; sub_x++)
                {
                    ash_x = ash_y + sub_x*subdiv_x;
                    na_x  = min(subdiv_x, na-(ash_x-ash));

                    fill_cell(nbs, grid, nbat,
                              ash_x, ash_x+na_x, atinfo, x,
                              grid->na_c*(cx*GPU_NSUBCELL_X+sub_x) + (dd_zone >> 2),
                              grid->na_c*(cy*GPU_NSUBCELL_Y+sub_y) + (dd_zone & 3),
                              grid->na_c*sub_z,
                              bb_work_aligned);
                }
            }
        }

        /* Set the unused atom indices to -1 */
        for (ind = na; ind < ncz*grid->na_sc; ind++)
        {
            nbs->a[ash+ind] = -1;
        }
    }
}

/* Determine in which grid column atoms should go */
static void calc_column_indices(nbnxn_grid_t *grid,
                                int a0, int a1,
                                rvec *x,
                                int dd_zone, const int *move,
                                int thread, int nthread,
                                int *cell,
                                int *cxy_na)
{
    int  n0, n1, i;
    int  cx, cy;

    /* We add one extra cell for particles which moved during DD */
    for (i = 0; i < grid->ncx*grid->ncy+1; i++)
    {
        cxy_na[i] = 0;
    }

    n0 = a0 + (int)((thread+0)*(a1 - a0))/nthread;
    n1 = a0 + (int)((thread+1)*(a1 - a0))/nthread;
    if (dd_zone == 0)
    {
        /* Home zone */
        for (i = n0; i < n1; i++)
        {
            if (move == NULL || move[i] >= 0)
            {
                /* We need to be careful with rounding,
                 * particles might be a few bits outside the local zone.
                 * The int cast takes care of the lower bound,
                 * we will explicitly take care of the upper bound.
                 */
                cx = (int)((x[i][XX] - grid->c0[XX])*grid->inv_sx);
                cy = (int)((x[i][YY] - grid->c0[YY])*grid->inv_sy);

#ifndef NDEBUG
                if (cx < 0 || cx > grid->ncx ||
                    cy < 0 || cy > grid->ncy)
                {
                    gmx_fatal(FARGS,
                              "grid cell cx %d cy %d out of range (max %d %d)\n"
                              "atom %f %f %f, grid->c0 %f %f",
                              cx, cy, grid->ncx, grid->ncy,
                              x[i][XX], x[i][YY], x[i][ZZ], grid->c0[XX], grid->c0[YY]);
                }
#endif
                /* Take care of potential rouding issues */
                cx = min(cx, grid->ncx - 1);
                cy = min(cy, grid->ncy - 1);

                /* For the moment cell will contain only the, grid local,
                 * x and y indices, not z.
                 */
                cell[i] = cx*grid->ncy + cy;
            }
            else
            {
                /* Put this moved particle after the end of the grid,
                 * so we can process it later without using conditionals.
                 */
                cell[i] = grid->ncx*grid->ncy;
            }

            cxy_na[cell[i]]++;
        }
    }
    else
    {
        /* Non-home zone */
        for (i = n0; i < n1; i++)
        {
            cx = (int)((x[i][XX] - grid->c0[XX])*grid->inv_sx);
            cy = (int)((x[i][YY] - grid->c0[YY])*grid->inv_sy);

            /* For non-home zones there could be particles outside
             * the non-bonded cut-off range, which have been communicated
             * for bonded interactions only. For the result it doesn't
             * matter where these end up on the grid. For performance
             * we put them in an extra row at the border.
             */
            cx = max(cx, 0);
            cx = min(cx, grid->ncx - 1);
            cy = max(cy, 0);
            cy = min(cy, grid->ncy - 1);

            /* For the moment cell will contain only the, grid local,
             * x and y indices, not z.
             */
            cell[i] = cx*grid->ncy + cy;

            cxy_na[cell[i]]++;
        }
    }
}

/* Determine in which grid cells the atoms should go */
static void calc_cell_indices(const nbnxn_search_t nbs,
                              int dd_zone,
                              nbnxn_grid_t *grid,
                              int a0, int a1,
                              const int *atinfo,
                              rvec *x,
                              const int *move,
                              nbnxn_atomdata_t *nbat)
{
    int   n0, n1, i;
    int   cx, cy, cxy, ncz_max = 0, ncz;
    int   nthread, thread;
    int  *cxy_na, cxy_na_i;

    nthread = gmx_omp_nthreads_get(emntPairsearch);

    /* For non-local zones we already have the indices */
    if (dd_zone == 0)
    {
#pragma omp parallel for num_threads(nthread) schedule(static)
        for (thread = 0; thread < nthread; thread++)
        {
            calc_column_indices(grid, a0, a1, x, dd_zone, move, thread, nthread,
                                nbs->cell, nbs->work[thread].cxy_na);
        }

        /* Make the cell index as a function of x and y */
        ncz_max          = 0;
        ncz              = 0;
        grid->cxy_ind[0] = 0;
        for (i = 0; i < grid->ncx*grid->ncy+1; i++)
        {
            /* We set ncz_max at the beginning of the loop iso at the end
             * to skip i=grid->ncx*grid->ncy which are moved particles
             * that do not need to be ordered on the grid.
             */
            if (ncz > ncz_max)
            {
                ncz_max = ncz;
            }
            cxy_na_i = nbs->work[0].cxy_na[i];
            for (thread = 1; thread < nthread; thread++)
            {
                cxy_na_i += nbs->work[thread].cxy_na[i];
            }
            ncz = (cxy_na_i + grid->na_sc - 1)/grid->na_sc;
            if (grid->na_cj == 2*grid->na_c)
            {
                /* Make the number of cell a multiple of 2 */
                ncz = (ncz + 1) & ~1;
            }
            grid->cxy_ind[i+1] = grid->cxy_ind[i] + ncz;
            /* Clear cxy_na, so we can reuse the array below */
            grid->cxy_na[i] = 0;
        }
        grid->nc = grid->cxy_ind[grid->ncx*grid->ncy] - grid->cxy_ind[0];
    }

    nbat->natoms = (grid->cell0 + grid->nc)*grid->na_sc;

    if (debug)
    {
        fprintf(debug, "ns na_sc %d na_c %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                grid->na_sc, grid->na_c, grid->nc,
                grid->ncx, grid->ncy, grid->nc/((double)(grid->ncx*grid->ncy)),
                dd_zone == 0 ? ncz_max : -1);
        if (gmx_debug_at)
        {
            i = 0;
            for (cy = 0; cy < grid->ncy; cy++)
            {
                for (cx = 0; cx < grid->ncx; cx++)
                {
                    fprintf(debug, " %2d", grid->cxy_ind[i+1]-grid->cxy_ind[i]);
                    i++;
                }
                fprintf(debug, "\n");
            }
        }
    }

    if (dd_zone == 0)
    {
        /* Make sure the work array for sorting is large enough */
        if (ncz_max*grid->na_sc*SGSF > nbs->work[0].sort_work_nalloc)
        {
            for (thread = 0; thread < nbs->nthread_max; thread++)
            {
                nbs->work[thread].sort_work_nalloc =
                    over_alloc_large(ncz_max*grid->na_sc*SGSF);
                srenew(nbs->work[thread].sort_work,
                       nbs->work[thread].sort_work_nalloc);
                /* When not in use, all elements should be -1 */
                for (i = 0; i < nbs->work[thread].sort_work_nalloc; i++)
                {
                    nbs->work[thread].sort_work[i] = -1;
                }
            }
        }

        /* Now we know the dimensions we can fill the grid.
         * This is the first, unsorted fill. We sort the columns after this.
         */
        for (i = a0; i < a1; i++)
        {
            /* At this point nbs->cell contains the local grid x,y indices */
            cxy = nbs->cell[i];
            nbs->a[(grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc + grid->cxy_na[cxy]++] = i;
        }

        /* Set the cell indices for the moved particles */
        n0 = grid->nc*grid->na_sc;
        n1 = grid->nc*grid->na_sc+grid->cxy_na[grid->ncx*grid->ncy];

        for (i = n0; i < n1; i++)
        {
            nbs->cell[nbs->a[i]] = i;
        }
    }
    else
    {
        /* Non-local zone: we have the exact cell indices already, use them */
        for (i = a0; i < a1; i++)
        {
            nbs->a[grid->cell0*grid->na_sc + nbs->cell[i]] = i;
        }
    }

    /* Sort the super-cell columns along z into the sub-cells.
     * For non-local zones there is only filling, no actual sorting.
     */
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        if (grid->bSimple)
        {
            sort_columns_simple(nbs, dd_zone, grid, a0, a1, atinfo, x, nbat,
                                ((thread+0)*grid->ncx*grid->ncy)/nthread,
                                ((thread+1)*grid->ncx*grid->ncy)/nthread,
                                nbs->work[thread].sort_work);
        }
        else
        {
            sort_columns_supersub(nbs, dd_zone, grid, a0, a1, atinfo, x, nbat,
                                  ((thread+0)*grid->ncx*grid->ncy)/nthread,
                                  ((thread+1)*grid->ncx*grid->ncy)/nthread,
                                  nbs->work[thread].sort_work);
        }
    }

    if (grid->bSimple && grid->na_cj == 2*grid->na_c)
    {
        combine_bounding_box_pairs(grid, grid->bb);
    }

    if (!grid->bSimple)
    {
        grid->nsubc_tot = 0;
        for (i = 0; i < grid->nc; i++)
        {
            grid->nsubc_tot += grid->nsubc[i];
        }
    }

    if (debug)
    {
        if (grid->bSimple)
        {
            print_bbsizes_simple(debug, grid);
        }
        else
        {
            fprintf(debug, "ns non-zero sub-cells: %d average atoms %.2f\n",
                    grid->nsubc_tot, (a1-a0)/(double)grid->nsubc_tot);

            print_bbsizes_supersub(debug, grid);
        }
    }
}

static void init_buffer_flags(nbnxn_buffer_flags_t *flags,
                              int                   natoms)
{
    int b;

    flags->nflag = (natoms + NBNXN_BUFFERFLAG_SIZE - 1)/NBNXN_BUFFERFLAG_SIZE;
    if (flags->nflag > flags->flag_nalloc)
    {
        flags->flag_nalloc = over_alloc_large(flags->nflag);
        srenew(flags->flag, flags->flag_nalloc);
    }
    for (b = 0; b < flags->nflag; b++)
    {
        bitmask_clear(&(flags->flag[b]));
    }
}

/* Set the atom order for a non-local grid for zone.
 * The grid sizes have been set up before in nbnxn_set_zone_grid.
 */
static void calc_grid_atom_order(nbnxn_search_t nbs, int zone,
                                 int at_start, int gmx_unused at_end)
{
    nbnxn_grid_t *grid;
    int           nat_tot_max, at, ci, cxy;

    grid = &nbs->grid[zone];

    nat_tot_max = (grid->cell0 + grid->nc)*grid->na_sc;

    if (nat_tot_max > nbs->cell_nalloc)
    {
        nbs->cell_nalloc = over_alloc_large(nat_tot_max);
        srenew(nbs->cell, nbs->cell_nalloc);
    }

    if (nat_tot_max > nbs->a_nalloc)
    {
        nbs->a_nalloc = over_alloc_large(nat_tot_max);
        srenew(nbs->a, nbs->a_nalloc);
    }

    at = at_start;
    ci = grid->cell0*grid->na_sc;
    for (cxy = 0; cxy < grid->ncx*grid->ncy; cxy++)
    {
        int na_column, i;

        /* Total number of atoms in the column, including fillers */
        na_column = (grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy])*grid->na_sc;

        ci = grid->cxy_ind[cxy]*grid->na_sc;
        for (i = 0; i < grid->cxy_na[cxy]; i++)
        {
            nbs->cell[at++] = ci++;
        }
    }

    assert(at == at_end);
}

/* Sets up a grid and puts the atoms on the grid.
 * This function only operates on one domain of the domain decompostion.
 * Note that without domain decomposition there is only one domain.
 */
void nbnxn_put_on_grid(nbnxn_search_t nbs,
                       int ePBC, matrix box,
                       int dd_zone,
                       rvec corner0, rvec corner1,
                       int a0, int a1,
                       real atom_density,
                       const int *atinfo,
                       rvec *x,
                       int nmoved, int *move,
                       int nb_kernel_type,
                       nbnxn_atomdata_t *nbat)
{
    nbnxn_grid_t *grid;
    int           n;
    int           nc_max_grid, nc_max;

    grid = &nbs->grid[dd_zone];

    nbs_cycle_start(&nbs->cc[enbsCCgrid]);

    grid->bSimple = nbnxn_kernel_pairlist_simple(nb_kernel_type);

    grid->na_c      = nbnxn_kernel_to_cluster_i_size(nb_kernel_type);
    grid->na_cj     = nbnxn_kernel_to_cluster_j_size(nb_kernel_type);
    grid->na_sc     = (grid->bSimple ? 1 : GPU_NSUBCELL)*grid->na_c;
    grid->na_c_2log = get_2log(grid->na_c);

    nbat->na_c = grid->na_c;

    if (dd_zone == 0)
    {
        grid->cell0 = 0;
    }
    else
    {
        grid->cell0 =
            (nbs->grid[dd_zone-1].cell0 + nbs->grid[dd_zone-1].nc)*
            nbs->grid[dd_zone-1].na_sc/grid->na_sc;
    }

    n = a1 - a0;

    if (dd_zone == 0)
    {
        nbs->ePBC = ePBC;
        copy_mat(box, nbs->box);

        /* Avoid zero density */
        if (atom_density > 0)
        {
            grid->atom_density = atom_density;
        }
        else
        {
            grid->atom_density = grid_atom_density(n-nmoved, corner0, corner1);
        }

        grid->cell0 = 0;

        nbs->natoms_local    = a1 - nmoved;
        /* We assume that nbnxn_put_on_grid is called first
         * for the local atoms (dd_zone=0).
         */
        nbs->natoms_nonlocal = a1 - nmoved;

        if (debug)
        {
            fprintf(debug, "natoms_local = %5d atom_density = %5.1f\n",
                    nbs->natoms_local, grid->atom_density);
        }
    }
    else
    {
        calc_grid_atom_order(nbs, dd_zone, a0, a1);

        nbs->natoms_nonlocal = max(nbs->natoms_nonlocal, a1);

        nc_max = grid->cell0 + grid->nc;
    }

    if (dd_zone == 0)
    {
        /* We always use the home zone (grid[0]) for setting the cell size,
         * since determining densities for non-local zones is difficult.
         */
        nc_max_grid = set_grid_size_xy(nbs, grid,
                                       n-nmoved, corner0, corner1,
                                       nbs->grid[0].atom_density);

        nc_max = grid->cell0 + nc_max_grid;

        if (a1 > nbs->cell_nalloc)
        {
            nbs->cell_nalloc = over_alloc_large(a1);
            srenew(nbs->cell, nbs->cell_nalloc);
        }

        /* To avoid conditionals we store the moved particles at the end of a,
         * make sure we have enough space.
         */
        if (nc_max*grid->na_sc + nmoved > nbs->a_nalloc)
        {
            nbs->a_nalloc = over_alloc_large(nc_max*grid->na_sc + nmoved);
            srenew(nbs->a, nbs->a_nalloc);
        }
    }

    grid_realloc_bb_flags(grid, nbs->bFEP, nc_max);

    /* We need padding up to a multiple of the buffer flag size: simply add */
    if (nc_max*grid->na_sc + NBNXN_BUFFERFLAG_SIZE > nbat->nalloc)
    {
        nbnxn_atomdata_realloc(nbat, nc_max*grid->na_sc+NBNXN_BUFFERFLAG_SIZE);
    }

    calc_cell_indices(nbs, dd_zone, grid, a0, a1, atinfo, x, move, nbat);

    if (dd_zone == 0)
    {
        nbat->natoms_local = nbat->natoms;
    }

    nbs_cycle_stop(&nbs->cc[enbsCCgrid]);
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nbnxn_search_t            nbs,
                                const gmx_domdec_zones_t *zones,
                                const int                *atinfo,
                                rvec                     *x,
                                int                       nb_kernel_type,
                                nbnxn_atomdata_t         *nbat)
{
    int  zone;
    rvec corner_unused;

    clear_rvec(corner_unused);

    for (zone = 1; zone < zones->n; zone++)
    {
        nbnxn_put_on_grid(nbs, nbs->ePBC, NULL,
                          zone, corner_unused, corner_unused,
                          zones->cg_range[zone],
                          zones->cg_range[zone+1],
                          -1,
                          atinfo,
                          x,
                          0, NULL,
                          nb_kernel_type,
                          nbat);
    }
}

/* Add simple grid type information to the local super/sub grid */
void nbnxn_grid_add_simple(nbnxn_search_t    nbs,
                           nbnxn_atomdata_t *nbat)
{
    nbnxn_grid_t *grid;
    float        *bbcz;
    nbnxn_bb_t   *bb;
    int           ncd, sc;
    int           nthreads gmx_unused;

    grid = &nbs->grid[0];

    if (grid->bSimple)
    {
        gmx_incons("nbnxn_grid_simple called with a simple grid");
    }

    ncd = grid->na_sc/NBNXN_CPU_CLUSTER_I_SIZE;

    if (grid->nc*ncd > grid->nc_nalloc_simple)
    {
        grid->nc_nalloc_simple = over_alloc_large(grid->nc*ncd);
        srenew(grid->bbcz_simple, grid->nc_nalloc_simple*NNBSBB_D);
        srenew(grid->bb_simple, grid->nc_nalloc_simple);
        srenew(grid->flags_simple, grid->nc_nalloc_simple);
        if (nbat->XFormat)
        {
            sfree_aligned(grid->bbj);
            snew_aligned(grid->bbj, grid->nc_nalloc_simple/2, 16);
        }
    }

    bbcz = grid->bbcz_simple;
    bb   = grid->bb_simple;

    nthreads = gmx_omp_nthreads_get(emntPairsearch);
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (sc = 0; sc < grid->nc; sc++)
    {
        int c, tx, na;

        for (c = 0; c < ncd; c++)
        {
            tx = sc*ncd + c;

            na = NBNXN_CPU_CLUSTER_I_SIZE;
            while (na > 0 &&
                   nbat->type[tx*NBNXN_CPU_CLUSTER_I_SIZE+na-1] == nbat->ntype-1)
            {
                na--;
            }

            if (na > 0)
            {
                switch (nbat->XFormat)
                {
                    case nbatX4:
                        /* PACK_X4==NBNXN_CPU_CLUSTER_I_SIZE, so this is simple */
                        calc_bounding_box_x_x4(na, nbat->x+tx*STRIDE_P4,
                                               bb+tx);
                        break;
                    case nbatX8:
                        /* PACK_X8>NBNXN_CPU_CLUSTER_I_SIZE, more complicated */
                        calc_bounding_box_x_x8(na, nbat->x+X8_IND_A(tx*NBNXN_CPU_CLUSTER_I_SIZE),
                                               bb+tx);
                        break;
                    default:
                        calc_bounding_box(na, nbat->xstride,
                                          nbat->x+tx*NBNXN_CPU_CLUSTER_I_SIZE*nbat->xstride,
                                          bb+tx);
                        break;
                }
                bbcz[tx*NNBSBB_D+0] = bb[tx].lower[BB_Z];
                bbcz[tx*NNBSBB_D+1] = bb[tx].upper[BB_Z];

                /* No interaction optimization yet here */
                grid->flags_simple[tx] = NBNXN_CI_DO_LJ(0) | NBNXN_CI_DO_COUL(0);
            }
            else
            {
                grid->flags_simple[tx] = 0;
            }
        }
    }

    if (grid->bSimple && grid->na_cj == 2*grid->na_c)
    {
        combine_bounding_box_pairs(grid, grid->bb_simple);
    }
}

void nbnxn_get_local_grid_sizes(nbnxn_search_t nbs,
                                int *ncx, int *ncy,
                                rvec *corner0, rvec *corner1,
                                real *column_size_x, real *column_size_y)
{
    const nbnxn_grid_t *grid;

    grid = &nbs->grid[0];

    *ncx = grid->ncx;
    *ncy = grid->ncy;
    if (corner0 != NULL)
    {
        copy_rvec(grid->c0, *corner0);
    }
    if (corner1 != NULL)
    {
        copy_rvec(grid->c1, *corner1);
    }
    if (column_size_x != NULL)
    {
        *column_size_x = grid->sx;
    }
    if (column_size_y != NULL)
    {
        *column_size_y = grid->sy;
    }
}

void nbnxn_get_local_grid_column(nbnxn_search_t nbs, int cx, int cy,
                                 nbnxn_bb_t *column_bb,
                                 int *bb_start,
                                 int *nbb, nbnxn_bb_t **bb, float **bbz,
                                 int *atom_start, int *bb_natoms, int *natoms)
{
    const nbnxn_grid_t *grid;
    int                 cxy;

    grid = &nbs->grid[0];

    assert(cx >= 0 && cx < grid->ncx);
    assert(cy >= 0 && cy < grid->ncy);

    column_bb->lower[XX] = grid->c0[XX] +  cx      * grid->sx;
    column_bb->upper[XX] = grid->c0[XX] + (cx + 1) * grid->sx;
    column_bb->lower[YY] = grid->c0[YY] +  cy      * grid->sy;
    column_bb->upper[YY] = grid->c0[YY] + (cy + 1) * grid->sy;
    column_bb->lower[ZZ] = grid->c0[ZZ];
    column_bb->upper[ZZ] = grid->c1[ZZ];

    cxy         = cx*grid->ncy + cy;

    if (grid->bSimple)
    {
        *bb_natoms = max(grid->na_c, grid->na_cj);

        *bb_start  = grid->cxy_ind[cxy];
        *nbb       = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        /* We return the largest of the x/y bounding boxes */
        if (grid->na_cj <= grid->na_c)
        {
            *bb    = grid->bb + grid->cxy_ind[cxy];
        }
        else
        {
            assert(grid->na_cj == 2*grid->na_c);
            /* j-clusters are twice as large as i, need to divide counts by 2 */
            *bb_start /= 2;
            *nbb      /= 2;
            *bb        = grid->bbj + grid->cxy_ind[cxy]/2;
        }
        *bbz   = NULL;
    }
    else
    {
        *bb_natoms = grid->na_sc;

        *bb_start  = grid->cxy_ind[cxy];
        *nbb       = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        *bb        = NULL;
        *bbz       = grid->bbcz + grid->cxy_ind[cxy]*NNBSBB_D;
    }

    *atom_start = nbs->a[grid->cxy_ind[cxy]*grid->na_sc];
    *natoms     = grid->cxy_na[cxy];
}

void nbnxn_set_zone_grid(nbnxn_search_t nbs,
                         int zone,
                         int ncx, int ncy,
                         rvec corner0, rvec corner1,
                         real column_size_x, real column_size_y,
                         const int *cxy_natoms)
{
    nbnxn_grid_t *grid;
    int           cxy;

    grid = &nbs->grid[zone];

    /* Copy the cluster sizes from the home zone */
    grid->na_c      = nbs->grid[0].na_c;
    grid->na_cj     = nbs->grid[0].na_cj;
    grid->na_sc     = nbs->grid[0].na_sc;
    grid->na_c_2log = nbs->grid[0].na_c_2log;

    /* The cell count continues after the previous grid/zone */
    grid->cell0 = nbs->grid[zone-1].cell0 + nbs->grid[zone-1].nc;

    grid->ncx     = ncx;
    grid->ncy     = ncy;
    copy_rvec(corner0, grid->c0);
    copy_rvec(corner1, grid->c1);
    rvec_sub(corner1, corner0, grid->size);
    grid->sx      = column_size_x;
    grid->sy      = column_size_y;
    grid->inv_sx  = 1/grid->sx;
    grid->inv_sy  = 1/grid->sy;

    if (grid->ncx*grid->ncy + 1 > grid->cxy_nalloc)
    {
        grid->cxy_nalloc = over_alloc_large(grid->ncx*grid->ncy + 1);
        srenew(grid->cxy_na, grid->cxy_nalloc);
        srenew(grid->cxy_ind, grid->cxy_nalloc);
    }

    grid->cxy_ind[0] = 0;
    for (cxy = 0; cxy < grid->ncx*grid->ncy; cxy++)
    {
        int nc_column;

        grid->cxy_na[cxy]    = cxy_natoms[cxy];
        nc_column            = (cxy_natoms[cxy] + grid->na_sc - 1)/grid->na_sc;
        if (grid->na_cj == 2*grid->na_c)
        {
            /* Make the number of cells a multiple of 2 */
            nc_column        = (nc_column + 1) & ~1;
        }
        grid->cxy_ind[cxy+1] = grid->cxy_ind[cxy] + nc_column;
    }
    grid->nc = grid->cxy_ind[grid->ncx*grid->ncy];
}

void nbnxn_get_atomorder(const nbnxn_search_t nbs, const int **a, int *n)
{
    const nbnxn_grid_t *grid;

    grid = &nbs->grid[0];

    /* Return the atom order for the home cell (index 0) */
    *a  = nbs->a;

    *n = grid->cxy_ind[grid->ncx*grid->ncy]*grid->na_sc;
}

void nbnxn_set_atomorder(nbnxn_search_t nbs)
{
    nbnxn_grid_t *grid;
    int           ao, cx, cy, cxy, cz, j;

    /* Set the atom order for the home cell (index 0) */
    grid = &nbs->grid[0];

    ao = 0;
    for (cx = 0; cx < grid->ncx; cx++)
    {
        for (cy = 0; cy < grid->ncy; cy++)
        {
            cxy = cx*grid->ncy + cy;
            j   = grid->cxy_ind[cxy]*grid->na_sc;
            for (cz = 0; cz < grid->cxy_na[cxy]; cz++)
            {
                nbs->a[j]     = ao;
                nbs->cell[ao] = j;
                ao++;
                j++;
            }
        }
    }
}
