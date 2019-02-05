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

#include "nbnxn_grid.h"

#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_util.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/mdtypes/forcerec.h" // only for GET_CGINFO_*
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

struct gmx_domdec_zones_t;

static real grid_atom_density(int        numAtoms,
                              const rvec lowerCorner,
                              const rvec upperCorner)
{
    rvec size;

    if (numAtoms == 0)
    {
        /* To avoid zero density we use a minimum of 1 atom */
        numAtoms = 1;
    }

    rvec_sub(upperCorner, lowerCorner, size);

    return numAtoms/(size[XX]*size[YY]*size[ZZ]);
}

static void set_grid_size_xy(const nbnxn_search   *nbs,
                             nbnxn_grid_t         *grid,
                             int                   ddZone,
                             int                   numAtoms,
                             const rvec            lowerCorner,
                             const rvec            upperCorner,
                             real                  atomDensity)
{
    rvec size;
    real tlen, tlen_x, tlen_y;

    rvec_sub(upperCorner, lowerCorner, size);

    if (numAtoms > grid->na_sc)
    {
        GMX_ASSERT(atomDensity > 0, "With one or more atoms, the density should be positive");

        /* target cell length */
        if (grid->bSimple)
        {
            /* To minimize the zero interactions, we should make
             * the largest of the i/j cell cubic.
             */
            int numAtomsInCell = std::max(grid->na_c, grid->na_cj);

            /* Approximately cubic cells */
            tlen   = std::cbrt(numAtomsInCell/atomDensity);
            tlen_x = tlen;
            tlen_y = tlen;
        }
        else
        {
            /* Approximately cubic sub cells */
            tlen   = std::cbrt(grid->na_c/atomDensity);
            tlen_x = tlen*c_gpuNumClusterPerCellX;
            tlen_y = tlen*c_gpuNumClusterPerCellY;
        }
        /* We round ncx and ncy down, because we get less cell pairs
         * in the nbsist when the fixed cell dimensions (x,y) are
         * larger than the variable one (z) than the other way around.
         */
        grid->numCells[XX] = std::max(1, static_cast<int>(size[XX]/tlen_x));
        grid->numCells[YY] = std::max(1, static_cast<int>(size[YY]/tlen_y));
    }
    else
    {
        grid->numCells[XX] = 1;
        grid->numCells[YY] = 1;
    }

    for (int d = 0; d < DIM - 1; d++)
    {
        grid->cellSize[d]    = size[d]/grid->numCells[d];
        grid->invCellSize[d] = 1/grid->cellSize[d];
    }

    if (ddZone > 0)
    {
        /* This is a non-home zone, add an extra row of cells
         * for particles communicated for bonded interactions.
         * These can be beyond the cut-off. It doesn't matter where
         * they end up on the grid, but for performance it's better
         * if they don't end up in cells that can be within cut-off range.
         */
        grid->numCells[XX]++;
        grid->numCells[YY]++;
    }

    /* We need one additional cell entry for particles moved by DD */
    int numCellsXY = grid->numCells[XX]*grid->numCells[YY];
    grid->cxy_na.resize(numCellsXY + 1);
    grid->cxy_ind.resize(numCellsXY + 2);

    for (nbnxn_search_work_t &work : nbs->work)
    {
        work.cxy_na.resize(numCellsXY + 1);
    }

    /* Worst case scenario of 1 atom in each last cell */
    int maxNumCells;
    if (grid->na_cj <= grid->na_c)
    {
        maxNumCells = numAtoms/grid->na_sc + numCellsXY;
    }
    else
    {
        maxNumCells = numAtoms/grid->na_sc + numCellsXY*grid->na_cj/grid->na_c;
    }

    grid->nsubc.resize(maxNumCells);
    grid->bbcz.resize(maxNumCells*NNBSBB_D);

    /* This resize also zeros the contents, this avoid possible
     * floating exceptions in SIMD with the unused bb elements.
     */
    if (grid->bSimple)
    {
        grid->bb.resize(maxNumCells);
    }
    else
    {
#if NBNXN_BBXXXX
        grid->pbb.resize(maxNumCells*c_gpuNumClusterPerCell/STRIDE_PBB*NNBSBB_XXXX);
#else
        grid->bb.resize(maxNumCells*c_gpuNumClusterPerCell);
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
            grid->bbjStorage.resize(maxNumCells*grid->na_c/grid->na_cj);
            grid->bbj = grid->bbjStorage;
        }
    }

    grid->flags.resize(maxNumCells);
    if (nbs->bFEP)
    {
        grid->fep.resize(maxNumCells*grid->na_sc/grid->na_c);
    }

    copy_rvec(lowerCorner, grid->c0);
    copy_rvec(upperCorner, grid->c1);
    copy_rvec(size,        grid->size);
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
                       bool gmx_unused relevantAtomsAreWithinGridBounds,
                       int *a, int n,
                       gmx::ArrayRef<const gmx::RVec> x,
                       real h0, real invh, int n_per_h,
                       gmx::ArrayRef<int> sort)
{
    if (n <= 1)
    {
        /* Nothing to do */
        return;
    }

    GMX_ASSERT(n <= n_per_h, "We require n <= n_per_h");

    /* Transform the inverse range height into the inverse hole height */
    invh *= n_per_h*SORT_GRID_OVERSIZE;

    /* Set nsort to the maximum possible number of holes used.
     * In worst case all n elements end up in the last bin.
     */
    int nsort = n_per_h*SORT_GRID_OVERSIZE + n;

    /* Determine the index range used, so we can limit it for the second pass */
    int zi_min = INT_MAX;
    int zi_max = -1;

    /* Sort the particles using a simple index sort */
    for (int i = 0; i < n; i++)
    {
        /* The cast takes care of float-point rounding effects below zero.
         * This code assumes particles are less than 1/SORT_GRID_OVERSIZE
         * times the box height out of the box.
         */
        int zi = static_cast<int>((x[a[i]][dim] - h0)*invh);

#ifndef NDEBUG
        /* As we can have rounding effect, we use > iso >= here */
        if (relevantAtomsAreWithinGridBounds &&
            (zi < 0 || (dd_zone == 0 && zi > n_per_h*SORT_GRID_OVERSIZE)))
        {
            gmx_fatal(FARGS, "(int)((x[%d][%c]=%f - %f)*%f) = %d, not in 0 - %d*%d\n",
                      a[i], 'x'+dim, x[a[i]][dim], h0, invh, zi,
                      n_per_h, SORT_GRID_OVERSIZE);
        }
#endif
        if (zi < 0)
        {
            zi = 0;
        }

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
            zi_min   = std::min(zi_min, zi);
            zi_max   = std::max(zi_max, zi);
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
                int cp  = sort[zi];
                int zim = zi + 1;
                while (sort[zim] >= 0)
                {
                    int tmp   = sort[zim];
                    sort[zim] = cp;
                    cp        = tmp;
                    zim++;
                }
                sort[zim] = cp;
                zi_max    = std::max(zi_max, zim);
            }
            sort[zi] = a[i];
            zi_max   = std::max(zi_max, zi);
        }
    }

    int c = 0;
    if (!Backwards)
    {
        for (int zi = 0; zi < nsort; zi++)
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
        for (int zi = zi_max; zi >= zi_min; zi--)
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

#if GMX_DOUBLE
#define R2F_D(x) ((float)((x) >= 0 ? ((1-GMX_FLOAT_EPS)*(x)) : ((1+GMX_FLOAT_EPS)*(x))))
#define R2F_U(x) ((float)((x) >= 0 ? ((1+GMX_FLOAT_EPS)*(x)) : ((1-GMX_FLOAT_EPS)*(x))))
#else
#define R2F_D(x) (x)
#define R2F_U(x) (x)
#endif

/* Coordinate order x,y,z, bb order xyz0 */
static void calc_bounding_box(int na, int stride, const real *x, nbnxn_bb_t *bb)
{
    int  i;
    real xl, xh, yl, yh, zl, zh;

    i  = 0;
    xl = x[i+XX];
    xh = x[i+XX];
    yl = x[i+YY];
    yh = x[i+YY];
    zl = x[i+ZZ];
    zh = x[i+ZZ];
    i += stride;
    for (int j = 1; j < na; j++)
    {
        xl = std::min(xl, x[i+XX]);
        xh = std::max(xh, x[i+XX]);
        yl = std::min(yl, x[i+YY]);
        yh = std::max(yh, x[i+YY]);
        zl = std::min(zl, x[i+ZZ]);
        zh = std::max(zh, x[i+ZZ]);
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
    real xl, xh, yl, yh, zl, zh;

    xl = x[XX*c_packX4];
    xh = x[XX*c_packX4];
    yl = x[YY*c_packX4];
    yh = x[YY*c_packX4];
    zl = x[ZZ*c_packX4];
    zh = x[ZZ*c_packX4];
    for (int j = 1; j < na; j++)
    {
        xl = std::min(xl, x[j+XX*c_packX4]);
        xh = std::max(xh, x[j+XX*c_packX4]);
        yl = std::min(yl, x[j+YY*c_packX4]);
        yh = std::max(yh, x[j+YY*c_packX4]);
        zl = std::min(zl, x[j+ZZ*c_packX4]);
        zh = std::max(zh, x[j+ZZ*c_packX4]);
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
    real xl, xh, yl, yh, zl, zh;

    xl = x[XX*c_packX8];
    xh = x[XX*c_packX8];
    yl = x[YY*c_packX8];
    yh = x[YY*c_packX8];
    zl = x[ZZ*c_packX8];
    zh = x[ZZ*c_packX8];
    for (int j = 1; j < na; j++)
    {
        xl = std::min(xl, x[j+XX*c_packX8]);
        xh = std::max(xh, x[j+XX*c_packX8]);
        yl = std::min(yl, x[j+YY*c_packX8]);
        yh = std::max(yh, x[j+YY*c_packX8]);
        zl = std::min(zl, x[j+ZZ*c_packX8]);
        zh = std::max(zh, x[j+ZZ*c_packX8]);
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
gmx_unused static void calc_bounding_box_x_x4_halves(int na, const real *x,
                                                     nbnxn_bb_t *bb, nbnxn_bb_t *bbj)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    calc_bounding_box_x_x4(std::min(na, 2), x, bbj);

    if (na > 2)
    {
        calc_bounding_box_x_x4(std::min(na-2, 2), x+(c_packX4>>1), bbj+1);
    }
    else
    {
        /* Set the "empty" bounding box to the same as the first one,
         * so we don't need to treat special cases in the rest of the code.
         */
#if NBNXN_SEARCH_BB_SIMD4
        store4(&bbj[1].lower[0], load4(&bbj[0].lower[0]));
        store4(&bbj[1].upper[0], load4(&bbj[0].upper[0]));
#else
        bbj[1] = bbj[0];
#endif
    }

#if NBNXN_SEARCH_BB_SIMD4
    store4(&bb->lower[0], min(load4(&bbj[0].lower[0]), load4(&bbj[1].lower[0])));
    store4(&bb->upper[0], max(load4(&bbj[0].upper[0]), load4(&bbj[1].upper[0])));
#else
    {
        int i;

        for (i = 0; i < NNBSBB_C; i++)
        {
            bb->lower[i] = std::min(bbj[0].lower[i], bbj[1].lower[i]);
            bb->upper[i] = std::max(bbj[0].upper[i], bbj[1].upper[i]);
        }
    }
#endif
}

#if NBNXN_SEARCH_BB_SIMD4

/* Coordinate order xyz, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx(int na, int stride, const real *x, float *bb)
{
    int  i;
    real xl, xh, yl, yh, zl, zh;

    i  = 0;
    xl = x[i+XX];
    xh = x[i+XX];
    yl = x[i+YY];
    yh = x[i+YY];
    zl = x[i+ZZ];
    zh = x[i+ZZ];
    i += stride;
    for (int j = 1; j < na; j++)
    {
        xl = std::min(xl, x[i+XX]);
        xh = std::max(xh, x[i+XX]);
        yl = std::min(yl, x[i+YY]);
        yh = std::max(yh, x[i+YY]);
        zl = std::min(zl, x[i+ZZ]);
        zh = std::max(zh, x[i+ZZ]);
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

#if NBNXN_SEARCH_SIMD4_FLOAT_X_BB

/* Coordinate order xyz?, bb order xyz0 */
static void calc_bounding_box_simd4(int na, const float *x, nbnxn_bb_t *bb)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    Simd4Float bb_0_S, bb_1_S;
    Simd4Float x_S;

    bb_0_S = load4(x);
    bb_1_S = bb_0_S;

    for (int i = 1; i < na; i++)
    {
        x_S    = load4(x+i*NNBSBB_C);
        bb_0_S = min(bb_0_S, x_S);
        bb_1_S = max(bb_1_S, x_S);
    }

    store4(&bb->lower[0], bb_0_S);
    store4(&bb->upper[0], bb_1_S);
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
static void combine_bounding_box_pairs(nbnxn_grid_t                    *grid,
                                       gmx::ArrayRef<const nbnxn_bb_t>  bb)
{
    // TODO: During SIMDv2 transition only some archs use namespace (remove when done)
    using namespace gmx;

    for (int i = 0; i < grid->numCells[XX]*grid->numCells[YY]; i++)
    {
        /* Starting bb in a column is expected to be 2-aligned */
        int sc2 = grid->cxy_ind[i]>>1;
        /* For odd numbers skip the last bb here */
        int nc2 = (grid->cxy_na[i]+3)>>(2+1);
        int c2;
        for (c2 = sc2; c2 < sc2+nc2; c2++)
        {
#if NBNXN_SEARCH_BB_SIMD4
            Simd4Float min_S, max_S;

            min_S = min(load4(&bb[c2*2+0].lower[0]),
                        load4(&bb[c2*2+1].lower[0]));
            max_S = max(load4(&bb[c2*2+0].upper[0]),
                        load4(&bb[c2*2+1].upper[0]));
            store4(&grid->bbj[c2].lower[0], min_S);
            store4(&grid->bbj[c2].upper[0], max_S);
#else
            for (int j = 0; j < NNBSBB_C; j++)
            {
                grid->bbj[c2].lower[j] = std::min(bb[c2*2+0].lower[j],
                                                  bb[c2*2+1].lower[j]);
                grid->bbj[c2].upper[j] = std::max(bb[c2*2+0].upper[j],
                                                  bb[c2*2+1].upper[j]);
            }
#endif
        }
        if (((grid->cxy_na[i]+3)>>2) & 1)
        {
            /* The bb count in this column is odd: duplicate the last bb */
            for (int j = 0; j < NNBSBB_C; j++)
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
    dvec ba;

    clear_dvec(ba);
    for (int c = 0; c < grid->nc; c++)
    {
        for (int d = 0; d < DIM; d++)
        {
            ba[d] += grid->bb[c].upper[d] - grid->bb[c].lower[d];
        }
    }
    dsvmul(1.0/grid->nc, ba, ba);

    real avgCellSizeZ = (grid->atom_density > 0 ? grid->na_c/(grid->atom_density*grid->cellSize[XX]*grid->cellSize[YY]) : 0.0);

    fprintf(fp, "ns bb: grid %4.2f %4.2f %4.2f abs %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            grid->cellSize[XX],
            grid->cellSize[YY],
            avgCellSizeZ,
            ba[XX], ba[YY], ba[ZZ],
            ba[XX]*grid->invCellSize[XX],
            ba[YY]*grid->invCellSize[YY],
            grid->atom_density > 0 ? ba[ZZ]/avgCellSizeZ : 0.0);
}

/* Prints the average bb size, used for debug output */
static void print_bbsizes_supersub(FILE                *fp,
                                   const nbnxn_grid_t  *grid)
{
    int  ns;
    dvec ba;

    clear_dvec(ba);
    ns = 0;
    for (int c = 0; c < grid->nc; c++)
    {
#if NBNXN_BBXXXX
        for (int s = 0; s < grid->nsubc[c]; s += STRIDE_PBB)
        {
            int cs_w = (c*c_gpuNumClusterPerCell + s)/STRIDE_PBB;
            for (int i = 0; i < STRIDE_PBB; i++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    ba[d] +=
                        grid->pbb[cs_w*NNBSBB_XXXX+(DIM+d)*STRIDE_PBB+i] -
                        grid->pbb[cs_w*NNBSBB_XXXX+     d *STRIDE_PBB+i];
                }
            }
        }
#else
        for (int s = 0; s < grid->nsubc[c]; s++)
        {
            int cs = c*c_gpuNumClusterPerCell + s;
            for (int d = 0; d < DIM; d++)
            {
                ba[d] += grid->bb[cs].upper[d] - grid->bb[cs].lower[d];
            }
        }
#endif
        ns += grid->nsubc[c];
    }
    dsvmul(1.0/ns, ba, ba);

    real avgClusterSizeZ = (grid->atom_density > 0 ? grid->na_sc/(grid->atom_density*grid->cellSize[XX]*grid->cellSize[YY]*c_gpuNumClusterPerCellZ) : 0.0);

    fprintf(fp, "ns bb: grid %4.2f %4.2f %4.2f abs %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            grid->cellSize[XX]/c_gpuNumClusterPerCellX,
            grid->cellSize[YY]/c_gpuNumClusterPerCellY,
            avgClusterSizeZ,
            ba[XX], ba[YY], ba[ZZ],
            ba[XX]*c_gpuNumClusterPerCellX*grid->invCellSize[XX],
            ba[YY]*c_gpuNumClusterPerCellY*grid->invCellSize[YY],
            grid->atom_density > 0 ? ba[ZZ]/avgClusterSizeZ : 0.0);
}

/* Set non-bonded interaction flags for the current cluster.
 * Sorts atoms on LJ coefficients: !=0 first, ==0 at the end.
 */
static void sort_cluster_on_flag(int                 numAtomsInCluster,
                                 int                 atomStart,
                                 int                 atomEnd,
                                 const int          *atinfo,
                                 gmx::ArrayRef<int>  order,
                                 int                *flags)
{
    constexpr int c_maxNumAtomsInCluster = 8;
    int           sort1[c_maxNumAtomsInCluster];
    int           sort2[c_maxNumAtomsInCluster];

    GMX_ASSERT(numAtomsInCluster <= c_maxNumAtomsInCluster, "Need to increase c_maxNumAtomsInCluster to support larger clusters");

    *flags = 0;

    int subc = 0;
    for (int s = atomStart; s < atomEnd; s += numAtomsInCluster)
    {
        /* Make lists for this (sub-)cell on atoms with and without LJ */
        int      n1         = 0;
        int      n2         = 0;
        gmx_bool haveQ      = FALSE;
        int      a_lj_max   = -1;
        for (int a = s; a < std::min(s + numAtomsInCluster, atomEnd); a++)
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

            if (2*n1 <= numAtomsInCluster)
            {
                /* Only sort when strictly necessary. Ordering particles
                 * Ordering particles can lead to less accurate summation
                 * due to rounding, both for LJ and Coulomb interactions.
                 */
                if (2*(a_lj_max - s) >= numAtomsInCluster)
                {
                    for (int i = 0; i < n1; i++)
                    {
                        order[atomStart + i]      = sort1[i];
                    }
                    for (int j = 0; j < n2; j++)
                    {
                        order[atomStart + n1 + j] = sort2[j];
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
static void fill_cell(nbnxn_search                  *nbs,
                      nbnxn_grid_t                  *grid,
                      nbnxn_atomdata_t              *nbat,
                      int                            atomStart,
                      int                            atomEnd,
                      const int                     *atinfo,
                      gmx::ArrayRef<const gmx::RVec> x,
                      nbnxn_bb_t gmx_unused         *bb_work_aligned)
{
    const int numAtoms = atomEnd - atomStart;

    if (grid->bSimple)
    {
        /* Note that non-local grids are already sorted.
         * Then sort_cluster_on_flag will only set the flags and the sorting
         * will not affect the atom order.
         */
        sort_cluster_on_flag(grid->na_c, atomStart, atomEnd, atinfo, nbs->a,
                             grid->flags.data() + (atomStart >> grid->na_c_2log) - grid->cell0);
    }

    if (nbs->bFEP)
    {
        /* Set the fep flag for perturbed atoms in this (sub-)cell */

        /* The grid-local cluster/(sub-)cell index */
        int cell = (atomStart >> grid->na_c_2log) - grid->cell0*(grid->bSimple ? 1 : c_gpuNumClusterPerCell);
        grid->fep[cell] = 0;
        for (int at = atomStart; at < atomEnd; at++)
        {
            if (nbs->a[at] >= 0 && GET_CGINFO_FEP(atinfo[nbs->a[at]]))
            {
                grid->fep[cell] |= (1 << (at - atomStart));
            }
        }
    }

    /* Now we have sorted the atoms, set the cell indices */
    for (int at = atomStart; at < atomEnd; at++)
    {
        nbs->cell[nbs->a[at]] = at;
    }

    copy_rvec_to_nbat_real(nbs->a.data() + atomStart, numAtoms, grid->na_c,
                           as_rvec_array(x.data()),
                           nbat->XFormat, nbat->x().data(), atomStart);

    if (nbat->XFormat == nbatX4)
    {
        /* Store the bounding boxes as xyz.xyz. */
        size_t      offset = (atomStart - grid->cell0*grid->na_sc) >> grid->na_c_2log;
        nbnxn_bb_t *bb_ptr = grid->bb.data() + offset;

#if GMX_SIMD && GMX_SIMD_REAL_WIDTH == 2
        if (2*grid->na_cj == grid->na_c)
        {
            calc_bounding_box_x_x4_halves(numAtoms, nbat->x().data() + atom_to_x_index<c_packX4>(atomStart), bb_ptr,
                                          grid->bbj.data() + offset*2);
        }
        else
#endif
        {
            calc_bounding_box_x_x4(numAtoms, nbat->x().data() + atom_to_x_index<c_packX4>(atomStart), bb_ptr);
        }
    }
    else if (nbat->XFormat == nbatX8)
    {
        /* Store the bounding boxes as xyz.xyz. */
        size_t      offset = (atomStart - grid->cell0*grid->na_sc) >> grid->na_c_2log;
        nbnxn_bb_t *bb_ptr = grid->bb.data() + offset;

        calc_bounding_box_x_x8(numAtoms, nbat->x().data() + atom_to_x_index<c_packX8>(atomStart), bb_ptr);
    }
#if NBNXN_BBXXXX
    else if (!grid->bSimple)
    {
        /* Store the bounding boxes in a format convenient
         * for SIMD4 calculations: xxxxyyyyzzzz...
         */
        float *pbb_ptr =
            grid->pbb.data() +
            ((atomStart - grid->cell0*grid->na_sc) >> (grid->na_c_2log + STRIDE_PBB_2LOG))*NNBSBB_XXXX +
            (((atomStart - grid->cell0*grid->na_sc) >> grid->na_c_2log) & (STRIDE_PBB - 1));

#if NBNXN_SEARCH_SIMD4_FLOAT_X_BB
        if (nbat->XFormat == nbatXYZQ)
        {
            calc_bounding_box_xxxx_simd4(numAtoms, nbat->x().data() + atomStart*nbat->xstride,
                                         bb_work_aligned, pbb_ptr);
        }
        else
#endif
        {
            calc_bounding_box_xxxx(numAtoms, nbat->xstride, nbat->x().data() + atomStart*nbat->xstride,
                                   pbb_ptr);
        }
        if (gmx_debug_at)
        {
            fprintf(debug, "cell %4d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    atomStart >> grid->na_c_2log,
                    pbb_ptr[0*STRIDE_PBB], pbb_ptr[3*STRIDE_PBB],
                    pbb_ptr[1*STRIDE_PBB], pbb_ptr[4*STRIDE_PBB],
                    pbb_ptr[2*STRIDE_PBB], pbb_ptr[5*STRIDE_PBB]);
        }
    }
#endif
    else
    {
        /* Store the bounding boxes as xyz.xyz. */
        nbnxn_bb_t *bb_ptr = grid->bb.data() + ((atomStart - grid->cell0*grid->na_sc) >> grid->na_c_2log);

        calc_bounding_box(numAtoms, nbat->xstride, nbat->x().data() + atomStart*nbat->xstride,
                          bb_ptr);

        if (gmx_debug_at)
        {
            int bbo = (atomStart - grid->cell0*grid->na_sc)/grid->na_c;
            fprintf(debug, "cell %4d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    atomStart >> grid->na_c_2log,
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
static void sort_columns_simple(nbnxn_search *nbs,
                                int dd_zone,
                                nbnxn_grid_t *grid,
                                int atomStart, int atomEnd,
                                const int *atinfo,
                                gmx::ArrayRef<const gmx::RVec> x,
                                nbnxn_atomdata_t *nbat,
                                int cxy_start, int cxy_end,
                                gmx::ArrayRef<int> sort_work)
{
    if (debug)
    {
        fprintf(debug, "cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0, cxy_start, cxy_end, atomStart, atomEnd);
    }

    const bool relevantAtomsAreWithinGridBounds = (grid->maxAtomGroupRadius == 0);

    /* Sort the atoms within each x,y column in 3 dimensions */
    for (int cxy = cxy_start; cxy < cxy_end; cxy++)
    {
        int na  = grid->cxy_na[cxy];
        int ncz = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        int ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ, FALSE, dd_zone,
                   relevantAtomsAreWithinGridBounds,
                   nbs->a.data() + ash, na, x,
                   grid->c0[ZZ],
                   1.0/grid->size[ZZ], ncz*grid->na_sc,
                   sort_work);

        /* Fill the ncz cells in this column */
        int cfilled = grid->cxy_ind[cxy];
        for (int cz = 0; cz < ncz; cz++)
        {
            int c     = grid->cxy_ind[cxy] + cz;

            int ash_c = ash + cz*grid->na_sc;
            int na_c  = std::min(grid->na_sc, na-(ash_c-ash));

            fill_cell(nbs, grid, nbat,
                      ash_c, ash_c+na_c, atinfo, x,
                      nullptr);

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
        for (int ind = na; ind < ncz*grid->na_sc; ind++)
        {
            nbs->a[ash+ind] = -1;
        }
    }
}

/* Spatially sort the atoms within one grid column */
static void sort_columns_supersub(nbnxn_search *nbs,
                                  int dd_zone,
                                  nbnxn_grid_t *grid,
                                  int atomStart, int atomEnd,
                                  const int *atinfo,
                                  gmx::ArrayRef<const gmx::RVec> x,
                                  nbnxn_atomdata_t *nbat,
                                  int cxy_start, int cxy_end,
                                  gmx::ArrayRef<int> sort_work)
{
    nbnxn_bb_t  bb_work_array[2];
    nbnxn_bb_t *bb_work_aligned = reinterpret_cast<nbnxn_bb_t *>((reinterpret_cast<std::size_t>(bb_work_array + 1)) & (~(static_cast<std::size_t>(15))));

    if (debug)
    {
        fprintf(debug, "cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0, cxy_start, cxy_end, atomStart, atomEnd);
    }

    const bool relevantAtomsAreWithinGridBounds = (grid->maxAtomGroupRadius == 0);

    int        subdiv_x = grid->na_c;
    int        subdiv_y = c_gpuNumClusterPerCellX*subdiv_x;
    int        subdiv_z = c_gpuNumClusterPerCellY*subdiv_y;

    /* Sort the atoms within each x,y column in 3 dimensions.
     * Loop over all columns on the x/y grid.
     */
    for (int cxy = cxy_start; cxy < cxy_end; cxy++)
    {
        int gridX            = cxy/grid->numCells[YY];
        int gridY            = cxy - gridX*grid->numCells[YY];

        int numAtomsInColumn = grid->cxy_na[cxy];
        int numCellsInColumn = grid->cxy_ind[cxy + 1] - grid->cxy_ind[cxy];
        int ash              = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ, FALSE, dd_zone,
                   relevantAtomsAreWithinGridBounds,
                   nbs->a.data() + ash, numAtomsInColumn, x,
                   grid->c0[ZZ],
                   1.0/grid->size[ZZ], numCellsInColumn*grid->na_sc,
                   sort_work);

        /* This loop goes over the cells and clusters along z at once */
        for (int sub_z = 0; sub_z < numCellsInColumn*c_gpuNumClusterPerCellZ; sub_z++)
        {
            int ash_z = ash + sub_z*subdiv_z;
            int na_z  = std::min(subdiv_z, numAtomsInColumn - (ash_z - ash));
            int cz    = -1;
            /* We have already sorted on z */

            if (sub_z % c_gpuNumClusterPerCellZ == 0)
            {
                cz = sub_z/c_gpuNumClusterPerCellZ;
                int cell = grid->cxy_ind[cxy] + cz;

                /* The number of atoms in this cell/super-cluster */
                int numAtoms = std::min(grid->na_sc, numAtomsInColumn - (ash_z - ash));

                grid->nsubc[cell] = std::min(c_gpuNumClusterPerCell,
                                             (numAtoms + grid->na_c - 1)/grid->na_c);

                /* Store the z-boundaries of the bounding box of the cell */
                grid->bbcz[cell*NNBSBB_D  ] = x[nbs->a[ash_z]][ZZ];
                grid->bbcz[cell*NNBSBB_D+1] = x[nbs->a[ash_z + numAtoms - 1]][ZZ];
            }

            if (c_gpuNumClusterPerCellY > 1)
            {
                /* Sort the atoms along y */
                sort_atoms(YY, (sub_z & 1) != 0, dd_zone,
                           relevantAtomsAreWithinGridBounds,
                           nbs->a.data() + ash_z, na_z, x,
                           grid->c0[YY] + gridY*grid->cellSize[YY],
                           grid->invCellSize[YY], subdiv_z,
                           sort_work);
            }

            for (int sub_y = 0; sub_y < c_gpuNumClusterPerCellY; sub_y++)
            {
                int ash_y = ash_z + sub_y*subdiv_y;
                int na_y  = std::min(subdiv_y, numAtomsInColumn - (ash_y - ash));

                if (c_gpuNumClusterPerCellX > 1)
                {
                    /* Sort the atoms along x */
                    sort_atoms(XX, ((cz*c_gpuNumClusterPerCellY + sub_y) & 1) != 0, dd_zone,
                               relevantAtomsAreWithinGridBounds,
                               nbs->a.data() + ash_y, na_y, x,
                               grid->c0[XX] + gridX*grid->cellSize[XX],
                               grid->invCellSize[XX], subdiv_y,
                               sort_work);
                }

                for (int sub_x = 0; sub_x < c_gpuNumClusterPerCellX; sub_x++)
                {
                    int ash_x = ash_y + sub_x*subdiv_x;
                    int na_x  = std::min(subdiv_x, numAtomsInColumn - (ash_x - ash));

                    fill_cell(nbs, grid, nbat,
                              ash_x, ash_x + na_x, atinfo, x,
                              bb_work_aligned);
                }
            }
        }

        /* Set the unused atom indices to -1 */
        for (int ind = numAtomsInColumn; ind < numCellsInColumn*grid->na_sc; ind++)
        {
            nbs->a[ash + ind] = -1;
        }
    }
}

/* Sets the cell index in the cell array for atom \p atomIndex and increments the atom count for the grid column */
static void setCellAndAtomCount(gmx::ArrayRef<int>  cell,
                                int                 cellIndex,
                                gmx::ArrayRef<int>  cxy_na,
                                int                 atomIndex)
{
    cell[atomIndex]    = cellIndex;
    cxy_na[cellIndex] += 1;
}

/* Determine in which grid column atoms should go */
static void calc_column_indices(nbnxn_grid_t *grid,
                                const gmx::UpdateGroupsCog *updateGroupsCog,
                                int atomStart, int atomEnd,
                                gmx::ArrayRef<const gmx::RVec> x,
                                int dd_zone, const int *move,
                                int thread, int nthread,
                                gmx::ArrayRef<int> cell,
                                gmx::ArrayRef<int> cxy_na)
{
    /* We add one extra cell for particles which moved during DD */
    for (int i = 0; i < grid->numCells[XX]*grid->numCells[YY] + 1; i++)
    {
        cxy_na[i] = 0;
    }

    int taskAtomStart = atomStart + static_cast<int>((thread + 0)*(atomEnd - atomStart))/nthread;
    int taskAtomEnd   = atomStart + static_cast<int>((thread + 1)*(atomEnd - atomStart))/nthread;
    if (dd_zone == 0)
    {
        /* Home zone */
        for (int i = taskAtomStart; i < taskAtomEnd; i++)
        {
            if (move == nullptr || move[i] >= 0)
            {

                const gmx::RVec &coord = (updateGroupsCog ? updateGroupsCog->cogForAtom(i) : x[i]);

                /* We need to be careful with rounding,
                 * particles might be a few bits outside the local zone.
                 * The int cast takes care of the lower bound,
                 * we will explicitly take care of the upper bound.
                 */
                int cx = static_cast<int>((coord[XX] - grid->c0[XX])*grid->invCellSize[XX]);
                int cy = static_cast<int>((coord[YY] - grid->c0[YY])*grid->invCellSize[YY]);

#ifndef NDEBUG
                if (cx < 0 || cx > grid->numCells[XX] ||
                    cy < 0 || cy > grid->numCells[YY])
                {
                    gmx_fatal(FARGS,
                              "grid cell cx %d cy %d out of range (max %d %d)\n"
                              "atom %f %f %f, grid->c0 %f %f",
                              cx, cy, grid->numCells[XX], grid->numCells[YY],
                              x[i][XX], x[i][YY], x[i][ZZ], grid->c0[XX], grid->c0[YY]);
                }
#endif
                /* Take care of potential rouding issues */
                cx = std::min(cx, grid->numCells[XX] - 1);
                cy = std::min(cy, grid->numCells[YY] - 1);

                /* For the moment cell will contain only the, grid local,
                 * x and y indices, not z.
                 */
                setCellAndAtomCount(cell, cx*grid->numCells[YY] + cy, cxy_na, i);
            }
            else
            {
                /* Put this moved particle after the end of the grid,
                 * so we can process it later without using conditionals.
                 */
                setCellAndAtomCount(cell, grid->numCells[XX]*grid->numCells[YY], cxy_na, i);
            }
        }
    }
    else
    {
        /* Non-home zone */
        for (int i = taskAtomStart; i < taskAtomEnd; i++)
        {
            int cx = static_cast<int>((x[i][XX] - grid->c0[XX])*grid->invCellSize[XX]);
            int cy = static_cast<int>((x[i][YY] - grid->c0[YY])*grid->invCellSize[YY]);

            /* For non-home zones there could be particles outside
             * the non-bonded cut-off range, which have been communicated
             * for bonded interactions only. For the result it doesn't
             * matter where these end up on the grid. For performance
             * we put them in an extra row at the border.
             */
            cx = std::max(cx, 0);
            cx = std::min(cx, grid->numCells[XX] - 1);
            cy = std::max(cy, 0);
            cy = std::min(cy, grid->numCells[YY] - 1);

            /* For the moment cell will contain only the, grid local,
             * x and y indices, not z.
             */
            setCellAndAtomCount(cell, cx*grid->numCells[YY] + cy, cxy_na, i);
        }
    }
}

/* Resizes grid and atom data which depend on the number of cells */
static void resizeForNumberOfCells(const nbnxn_grid_t &grid,
                                   int                 numAtomsMoved,
                                   nbnxn_search       *nbs,
                                   nbnxn_atomdata_t   *nbat)
{
    int numNbnxnAtoms = (grid.cell0 + grid.nc)*grid.na_sc;

    /* Note: nbs->cell was already resized before */

    /* To avoid conditionals we store the moved particles at the end of a,
     * make sure we have enough space.
     */
    nbs->a.resize(numNbnxnAtoms + numAtomsMoved);

    /* Make space in nbat for storing the atom coordinates */
    nbat->resizeCoordinateBuffer(numNbnxnAtoms);
}

/* Determine in which grid cells the atoms should go */
static void
calc_cell_indices(nbnxn_search                   *nbs,
                  int                             ddZone,
                  nbnxn_grid_t                   *grid,
                  const gmx::UpdateGroupsCog     *updateGroupsCog,
                  int                             atomStart,
                  int                             atomEnd,
                  const int                      *atinfo,
                  gmx::ArrayRef<const gmx::RVec>  x,
                  int                             numAtomsMoved,
                  const int                      *move,
                  nbnxn_atomdata_t               *nbat)
{
    /* First compute all grid/column indices and store them in nbs->cell */
    nbs->cell.resize(atomEnd);

    const int nthread = gmx_omp_nthreads_get(emntPairsearch);

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            calc_column_indices(grid, updateGroupsCog,
                                atomStart, atomEnd, x,
                                ddZone, move, thread, nthread,
                                nbs->cell, nbs->work[thread].cxy_na);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    /* Make the cell index as a function of x and y */
    int ncz_max      = 0;
    int ncz          = 0;
    grid->cxy_ind[0] = 0;
    for (int i = 0; i < grid->numCells[XX]*grid->numCells[YY] + 1; i++)
    {
        /* We set ncz_max at the beginning of the loop iso at the end
         * to skip i=grid->ncx*grid->numCells[YY] which are moved particles
         * that do not need to be ordered on the grid.
         */
        if (ncz > ncz_max)
        {
            ncz_max = ncz;
        }
        int cxy_na_i = nbs->work[0].cxy_na[i];
        for (int thread = 1; thread < nthread; thread++)
        {
            cxy_na_i += nbs->work[thread].cxy_na[i];
        }
        ncz = (cxy_na_i + grid->na_sc - 1)/grid->na_sc;
        if (nbat->XFormat == nbatX8)
        {
            /* Make the number of cell a multiple of 2 */
            ncz = (ncz + 1) & ~1;
        }
        grid->cxy_ind[i+1] = grid->cxy_ind[i] + ncz;
        /* Clear cxy_na, so we can reuse the array below */
        grid->cxy_na[i] = 0;
    }
    grid->nc = grid->cxy_ind[grid->numCells[XX]*grid->numCells[YY]] - grid->cxy_ind[0];

    resizeForNumberOfCells(*grid, numAtomsMoved, nbs, nbat);

    if (debug)
    {
        fprintf(debug, "ns na_sc %d na_c %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                grid->na_sc, grid->na_c, grid->nc,
                grid->numCells[XX], grid->numCells[YY], grid->nc/(static_cast<double>(grid->numCells[XX]*grid->numCells[YY])),
                ncz_max);
        if (gmx_debug_at)
        {
            int i = 0;
            for (int cy = 0; cy < grid->numCells[YY]; cy++)
            {
                for (int cx = 0; cx < grid->numCells[XX]; cx++)
                {
                    fprintf(debug, " %2d", grid->cxy_ind[i+1]-grid->cxy_ind[i]);
                    i++;
                }
                fprintf(debug, "\n");
            }
        }
    }

    /* Make sure the work array for sorting is large enough */
    if (ncz_max*grid->na_sc*SGSF > gmx::index(nbs->work[0].sortBuffer.size()))
    {
        for (nbnxn_search_work_t &work : nbs->work)
        {
            /* Elements not in use should be -1 */
            work.sortBuffer.resize(ncz_max*grid->na_sc*SGSF, -1);
        }
    }

    /* Now we know the dimensions we can fill the grid.
     * This is the first, unsorted fill. We sort the columns after this.
     */
    for (int i = atomStart; i < atomEnd; i++)
    {
        /* At this point nbs->cell contains the local grid x,y indices */
        int cxy = nbs->cell[i];
        nbs->a[(grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc + grid->cxy_na[cxy]++] = i;
    }

    if (ddZone == 0)
    {
        /* Set the cell indices for the moved particles */
        int n0 = grid->nc*grid->na_sc;
        int n1 = grid->nc*grid->na_sc+grid->cxy_na[grid->numCells[XX]*grid->numCells[YY]];
        if (ddZone == 0)
        {
            for (int i = n0; i < n1; i++)
            {
                nbs->cell[nbs->a[i]] = i;
            }
        }
    }

    /* Sort the super-cell columns along z into the sub-cells. */
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            int columnStart = ((thread + 0)*grid->numCells[XX]*grid->numCells[YY])/nthread;
            int columnEnd   = ((thread + 1)*grid->numCells[XX]*grid->numCells[YY])/nthread;
            if (grid->bSimple)
            {
                sort_columns_simple(nbs, ddZone, grid, atomStart, atomEnd, atinfo, x, nbat,
                                    columnStart, columnEnd,
                                    nbs->work[thread].sortBuffer);
            }
            else
            {
                sort_columns_supersub(nbs, ddZone, grid, atomStart, atomEnd, atinfo, x, nbat,
                                      columnStart, columnEnd,
                                      nbs->work[thread].sortBuffer);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    if (grid->bSimple && nbat->XFormat == nbatX8)
    {
        combine_bounding_box_pairs(grid, grid->bb);
    }

    if (!grid->bSimple)
    {
        grid->nsubc_tot = 0;
        for (int i = 0; i < grid->nc; i++)
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
                    grid->nsubc_tot, (atomEnd - atomStart)/static_cast<double>(grid->nsubc_tot));

            print_bbsizes_supersub(debug, grid);
        }
    }
}

/* Sets up a grid and puts the atoms on the grid.
 * This function only operates on one domain of the domain decompostion.
 * Note that without domain decomposition there is only one domain.
 */
void nbnxn_put_on_grid(nbnxn_search_t                  nbs,
                       int                             ePBC,
                       const matrix                    box,
                       int                             ddZone,
                       const rvec                      lowerCorner,
                       const rvec                      upperCorner,
                       const gmx::UpdateGroupsCog     *updateGroupsCog,
                       int                             atomStart,
                       int                             atomEnd,
                       real                            atomDensity,
                       const int                      *atinfo,
                       gmx::ArrayRef<const gmx::RVec>  x,
                       int                             numAtomsMoved,
                       const int                      *move,
                       int                             nb_kernel_type,
                       nbnxn_atomdata_t               *nbat)
{
    nbnxn_grid_t *grid = &nbs->grid[ddZone];

    nbs_cycle_start(&nbs->cc[enbsCCgrid]);

    grid->bSimple = nbnxn_kernel_pairlist_simple(nb_kernel_type);

    grid->na_c      = nbnxn_kernel_to_cluster_i_size(nb_kernel_type);
    grid->na_cj     = nbnxn_kernel_to_cluster_j_size(nb_kernel_type);
    grid->na_sc     = (grid->bSimple ? 1 : c_gpuNumClusterPerCell)*grid->na_c;
    grid->na_c_2log = get_2log(grid->na_c);

    if (ddZone == 0)
    {
        grid->cell0 = 0;
    }
    else
    {
        grid->cell0 =
            (nbs->grid[ddZone - 1].cell0 + nbs->grid[ddZone - 1].nc)*
            nbs->grid[ddZone- 1].na_sc/grid->na_sc;
    }

    const int n = atomEnd - atomStart;

    if (ddZone == 0)
    {
        nbs->ePBC = ePBC;
        copy_mat(box, nbs->box);

        /* Avoid zero density */
        if (atomDensity > 0)
        {
            grid->atom_density = atomDensity;
        }
        else
        {
            grid->atom_density = grid_atom_density(n - numAtomsMoved, lowerCorner, upperCorner);
        }

        grid->cell0 = 0;

        nbs->natoms_local    = atomEnd - numAtomsMoved;
        /* We assume that nbnxn_put_on_grid is called first
         * for the local atoms (ddZone=0).
         */
        nbs->natoms_nonlocal = atomEnd - numAtomsMoved;

        /* When not using atom groups, all atoms should be within the grid bounds */
        grid->maxAtomGroupRadius = (updateGroupsCog ? updateGroupsCog->maxUpdateGroupRadius() : 0);
        /* For the non-local grids the situation is the same as for the local */
        for (size_t g = 1; g < nbs->grid.size(); g++)
        {
            nbs->grid[g].maxAtomGroupRadius = grid->maxAtomGroupRadius;
        }

        if (debug)
        {
            fprintf(debug, "natoms_local = %5d atom_density = %5.1f\n",
                    nbs->natoms_local, grid->atom_density);
        }
    }
    else
    {
        nbs->natoms_nonlocal = std::max(nbs->natoms_nonlocal, atomEnd);
    }

    /* We always use the home zone (grid[0]) for setting the cell size,
     * since determining densities for non-local zones is difficult.
     */
    set_grid_size_xy(nbs, grid,
                     ddZone, n - numAtomsMoved,
                     lowerCorner, upperCorner,
                     nbs->grid[0].atom_density);

    calc_cell_indices(nbs, ddZone, grid, updateGroupsCog, atomStart, atomEnd, atinfo, x, numAtomsMoved, move, nbat);

    if (ddZone == 0)
    {
        nbat->natoms_local = nbat->numAtoms();
    }
    if (ddZone == gmx::ssize(nbs->grid) - 1)
    {
        /* We are done setting up all grids, we can resize the force buffers */
        nbat->resizeForceBuffers();
    }

    nbs_cycle_stop(&nbs->cc[enbsCCgrid]);
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nbnxn_search_t                   nbs,
                                const struct gmx_domdec_zones_t *zones,
                                const int                       *atinfo,
                                gmx::ArrayRef<const gmx::RVec>   x,
                                int                              nb_kernel_type,
                                nbnxn_atomdata_t                *nbat)
{
    for (int zone = 1; zone < zones->n; zone++)
    {
        rvec c0, c1;
        for (int d = 0; d < DIM; d++)
        {
            c0[d] = zones->size[zone].bb_x0[d];
            c1[d] = zones->size[zone].bb_x1[d];
        }

        nbnxn_put_on_grid(nbs, nbs->ePBC, nullptr,
                          zone, c0, c1,
                          nullptr,
                          zones->cg_range[zone],
                          zones->cg_range[zone+1],
                          -1,
                          atinfo,
                          x,
                          0, nullptr,
                          nb_kernel_type,
                          nbat);
    }
}

void nbnxn_get_ncells(nbnxn_search_t nbs, int *ncx, int *ncy)
{
    *ncx = nbs->grid[0].numCells[XX];
    *ncy = nbs->grid[0].numCells[YY];
}

gmx::ArrayRef<const int> nbnxn_get_atomorder(const nbnxn_search *nbs)
{
    /* Return the atom order for the home cell (index 0) */
    const nbnxn_grid_t &grid       = nbs->grid[0];

    int                 numIndices = grid.cxy_ind[grid.numCells[XX]*grid.numCells[YY]]*grid.na_sc;

    return gmx::constArrayRefFromArray(nbs->a.data(), numIndices);
}

void nbnxn_set_atomorder(nbnxn_search_t nbs)
{
    /* Set the atom order for the home cell (index 0) */
    nbnxn_grid_t *grid = &nbs->grid[0];

    int           ao = 0;
    for (int cx = 0; cx < grid->numCells[XX]; cx++)
    {
        for (int cy = 0; cy < grid->numCells[YY]; cy++)
        {
            int cxy = cx*grid->numCells[YY] + cy;
            int j   = grid->cxy_ind[cxy]*grid->na_sc;
            for (int cz = 0; cz < grid->cxy_na[cxy]; cz++)
            {
                nbs->a[j]     = ao;
                nbs->cell[ao] = j;
                ao++;
                j++;
            }
        }
    }
}

gmx::ArrayRef<const int> nbnxn_get_gridindices(const nbnxn_search *nbs)
{
    return nbs->cell;
}
