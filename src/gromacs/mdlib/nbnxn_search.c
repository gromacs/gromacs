/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "nbnxn_search.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/smalloc.h"

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

#ifdef GMX_NBNXN_SIMD

/* The functions below are macros as they are performance sensitive */

/* 4x4 list, pack=4: no complex conversion required */
/* i-cluster to j-cluster conversion */
#define CI_TO_CJ_J4(ci)   (ci)
/* cluster index to coordinate array index conversion */
#define X_IND_CI_J4(ci)  ((ci)*STRIDE_P4)
#define X_IND_CJ_J4(cj)  ((cj)*STRIDE_P4)

/* 4x2 list, pack=4: j-cluster size is half the packing width */
/* i-cluster to j-cluster conversion */
#define CI_TO_CJ_J2(ci)  ((ci)<<1)
/* cluster index to coordinate array index conversion */
#define X_IND_CI_J2(ci)  ((ci)*STRIDE_P4)
#define X_IND_CJ_J2(cj)  (((cj)>>1)*STRIDE_P4 + ((cj) & 1)*(PACK_X4>>1))

/* 4x8 list, pack=8: i-cluster size is half the packing width */
/* i-cluster to j-cluster conversion */
#define CI_TO_CJ_J8(ci)  ((ci)>>1)
/* cluster index to coordinate array index conversion */
#define X_IND_CI_J8(ci)  (((ci)>>1)*STRIDE_P8 + ((ci) & 1)*(PACK_X8>>1))
#define X_IND_CJ_J8(cj)  ((cj)*STRIDE_P8)

/* The j-cluster size is matched to the SIMD width */
#if GMX_SIMD_REAL_WIDTH == 2
#define CI_TO_CJ_SIMD_4XN(ci)  CI_TO_CJ_J2(ci)
#define X_IND_CI_SIMD_4XN(ci)  X_IND_CI_J2(ci)
#define X_IND_CJ_SIMD_4XN(cj)  X_IND_CJ_J2(cj)
#else
#if GMX_SIMD_REAL_WIDTH == 4
#define CI_TO_CJ_SIMD_4XN(ci)  CI_TO_CJ_J4(ci)
#define X_IND_CI_SIMD_4XN(ci)  X_IND_CI_J4(ci)
#define X_IND_CJ_SIMD_4XN(cj)  X_IND_CJ_J4(cj)
#else
#if GMX_SIMD_REAL_WIDTH == 8
#define CI_TO_CJ_SIMD_4XN(ci)  CI_TO_CJ_J8(ci)
#define X_IND_CI_SIMD_4XN(ci)  X_IND_CI_J8(ci)
#define X_IND_CJ_SIMD_4XN(cj)  X_IND_CJ_J8(cj)
/* Half SIMD with j-cluster size */
#define CI_TO_CJ_SIMD_2XNN(ci) CI_TO_CJ_J4(ci)
#define X_IND_CI_SIMD_2XNN(ci) X_IND_CI_J4(ci)
#define X_IND_CJ_SIMD_2XNN(cj) X_IND_CJ_J4(cj)
#else
#if GMX_SIMD_REAL_WIDTH == 16
#define CI_TO_CJ_SIMD_2XNN(ci) CI_TO_CJ_J8(ci)
#define X_IND_CI_SIMD_2XNN(ci) X_IND_CI_J8(ci)
#define X_IND_CJ_SIMD_2XNN(cj) X_IND_CJ_J8(cj)
#else
#error "unsupported GMX_SIMD_REAL_WIDTH"
#endif
#endif
#endif
#endif

#endif /* GMX_NBNXN_SIMD */


#ifdef NBNXN_SEARCH_BB_SIMD4
/* Store bounding boxes corners as quadruplets: xxxxyyyyzzzz */
#define NBNXN_BBXXXX
/* Size of bounding box corners quadruplet */
#define NNBSBB_XXXX      (NNBSBB_D*DIM*STRIDE_PBB)
#endif

/* We shift the i-particles backward for PBC.
 * This leads to more conditionals than shifting forward.
 * We do this to get more balanced pair lists.
 */
#define NBNXN_SHIFT_BACKWARD


/* This define is a lazy way to avoid interdependence of the grid
 * and searching data structures.
 */
#define NBNXN_NA_SC_MAX (GPU_NSUBCELL*NBNXN_GPU_CLUSTER_SIZE)


static void nbs_cycle_clear(nbnxn_cycle_t *cc)
{
    int i;

    for (i = 0; i < enbsCCnr; i++)
    {
        cc[i].count = 0;
        cc[i].c     = 0;
    }
}

static double Mcyc_av(const nbnxn_cycle_t *cc)
{
    return (double)cc->c*1e-6/cc->count;
}

static void nbs_cycle_print(FILE *fp, const nbnxn_search_t nbs)
{
    int n;
    int t;

    fprintf(fp, "\n");
    fprintf(fp, "ns %4d grid %4.1f search %4.1f red.f %5.3f",
            nbs->cc[enbsCCgrid].count,
            Mcyc_av(&nbs->cc[enbsCCgrid]),
            Mcyc_av(&nbs->cc[enbsCCsearch]),
            Mcyc_av(&nbs->cc[enbsCCreducef]));

    if (nbs->nthread_max > 1)
    {
        if (nbs->cc[enbsCCcombine].count > 0)
        {
            fprintf(fp, " comb %5.2f",
                    Mcyc_av(&nbs->cc[enbsCCcombine]));
        }
        fprintf(fp, " s. th");
        for (t = 0; t < nbs->nthread_max; t++)
        {
            fprintf(fp, " %4.1f",
                    Mcyc_av(&nbs->work[t].cc[enbsCCsearch]));
        }
    }
    fprintf(fp, "\n");
}

static void nbnxn_grid_init(nbnxn_grid_t * grid)
{
    grid->cxy_na      = NULL;
    grid->cxy_ind     = NULL;
    grid->cxy_nalloc  = 0;
    grid->bb          = NULL;
    grid->bbj         = NULL;
    grid->nc_nalloc   = 0;
}

static int get_2log(int n)
{
    int log2;

    log2 = 0;
    while ((1<<log2) < n)
    {
        log2++;
    }
    if ((1<<log2) != n)
    {
        gmx_fatal(FARGS, "nbnxn na_c (%d) is not a power of 2", n);
    }

    return log2;
}

int nbnxn_kernel_to_ci_size(int nb_kernel_type)
{
    switch (nb_kernel_type)
    {
        case nbnxnk4x4_PlainC:
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
            return NBNXN_CPU_CLUSTER_I_SIZE;
        case nbnxnk8x8x8_GPU:
        case nbnxnk8x8x8_PlainC:
            /* The cluster size for super/sub lists is only set here.
             * Any value should work for the pair-search and atomdata code.
             * The kernels, of course, might require a particular value.
             */
            return NBNXN_GPU_CLUSTER_SIZE;
        default:
            gmx_incons("unknown kernel type");
    }

    return 0;
}

int nbnxn_kernel_to_cj_size(int nb_kernel_type)
{
    int nbnxn_simd_width = 0;
    int cj_size          = 0;

#ifdef GMX_NBNXN_SIMD
    nbnxn_simd_width = GMX_SIMD_REAL_WIDTH;
#endif

    switch (nb_kernel_type)
    {
        case nbnxnk4x4_PlainC:
            cj_size = NBNXN_CPU_CLUSTER_I_SIZE;
            break;
        case nbnxnk4xN_SIMD_4xN:
            cj_size = nbnxn_simd_width;
            break;
        case nbnxnk4xN_SIMD_2xNN:
            cj_size = nbnxn_simd_width/2;
            break;
        case nbnxnk8x8x8_GPU:
        case nbnxnk8x8x8_PlainC:
            cj_size = nbnxn_kernel_to_ci_size(nb_kernel_type);
            break;
        default:
            gmx_incons("unknown kernel type");
    }

    return cj_size;
}

static int ci_to_cj(int na_cj_2log, int ci)
{
    switch (na_cj_2log)
    {
        case 2: return ci;     break;
        case 1: return (ci<<1); break;
        case 3: return (ci>>1); break;
    }

    return 0;
}

gmx_bool nbnxn_kernel_pairlist_simple(int nb_kernel_type)
{
    if (nb_kernel_type == nbnxnkNotSet)
    {
        gmx_fatal(FARGS, "Non-bonded kernel type not set for Verlet-style pair-list.");
    }

    switch (nb_kernel_type)
    {
        case nbnxnk8x8x8_GPU:
        case nbnxnk8x8x8_PlainC:
            return FALSE;

        case nbnxnk4x4_PlainC:
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
            return TRUE;

        default:
            gmx_incons("Invalid nonbonded kernel type passed!");
            return FALSE;
    }
}

/* Initializes a single nbnxn_pairlist_t data structure */
static void nbnxn_init_pairlist_fep(t_nblist *nl)
{
    nl->type        = GMX_NBLIST_INTERACTION_FREE_ENERGY;
    nl->igeometry   = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
    /* The interaction functions are set in the free energy kernel fuction */
    nl->ivdw        = -1;
    nl->ivdwmod     = -1;
    nl->ielec       = -1;
    nl->ielecmod    = -1;

    nl->maxnri      = 0;
    nl->maxnrj      = 0;
    nl->nri         = 0;
    nl->nrj         = 0;
    nl->iinr        = NULL;
    nl->gid         = NULL;
    nl->shift       = NULL;
    nl->jindex      = NULL;
    nl->jjnr        = NULL;
    nl->excl_fep    = NULL;

}

void nbnxn_init_search(nbnxn_search_t    * nbs_ptr,
                       ivec               *n_dd_cells,
                       gmx_domdec_zones_t *zones,
                       gmx_bool            bFEP,
                       int                 nthread_max)
{
    nbnxn_search_t nbs;
    int            d, g, t;

    snew(nbs, 1);
    *nbs_ptr = nbs;

    nbs->bFEP   = bFEP;

    nbs->DomDec = (n_dd_cells != NULL);

    clear_ivec(nbs->dd_dim);
    nbs->ngrid = 1;
    if (nbs->DomDec)
    {
        nbs->zones = zones;

        for (d = 0; d < DIM; d++)
        {
            if ((*n_dd_cells)[d] > 1)
            {
                nbs->dd_dim[d] = 1;
                /* Each grid matches a DD zone */
                nbs->ngrid *= 2;
            }
        }
    }

    snew(nbs->grid, nbs->ngrid);
    for (g = 0; g < nbs->ngrid; g++)
    {
        nbnxn_grid_init(&nbs->grid[g]);
    }
    nbs->cell        = NULL;
    nbs->cell_nalloc = 0;
    nbs->a           = NULL;
    nbs->a_nalloc    = 0;

    nbs->nthread_max = nthread_max;

    /* Initialize the work data structures for each thread */
    snew(nbs->work, nbs->nthread_max);
    for (t = 0; t < nbs->nthread_max; t++)
    {
        nbs->work[t].cxy_na           = NULL;
        nbs->work[t].cxy_na_nalloc    = 0;
        nbs->work[t].sort_work        = NULL;
        nbs->work[t].sort_work_nalloc = 0;

        snew(nbs->work[t].nbl_fep, 1);
        nbnxn_init_pairlist_fep(nbs->work[t].nbl_fep);
    }

    /* Initialize detailed nbsearch cycle counting */
    nbs->print_cycles = (getenv("GMX_NBNXN_CYCLE") != 0);
    nbs->search_count = 0;
    nbs_cycle_clear(nbs->cc);
    for (t = 0; t < nbs->nthread_max; t++)
    {
        nbs_cycle_clear(nbs->work[t].cc);
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

static int set_grid_size_xy(const nbnxn_search_t nbs,
                            nbnxn_grid_t *grid,
                            int dd_zone,
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
        /* We round ncx and ncy down, because we get less cell pairs
         * in the nbsist when the fixed cell dimensions (x,y) are
         * larger than the variable one (z) than the other way around.
         */
        grid->ncx = max(1, (int)(size[XX]/tlen_x));
        grid->ncy = max(1, (int)(size[YY]/tlen_y));
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

    if (dd_zone > 0)
    {
        /* This is a non-home zone, add an extra row of cells
         * for particles communicated for bonded interactions.
         * These can be beyond the cut-off. It doesn't matter where
         * they end up on the grid, but for performance it's better
         * if they don't end up in cells that can be within cut-off range.
         */
        grid->ncx++;
        grid->ncy++;
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
        if (nbs->bFEP)
        {
            srenew(grid->fep, grid->nc_nalloc*grid->na_sc/grid->na_c);
        }
    }

    copy_rvec(corner0, grid->c0);
    copy_rvec(corner1, grid->c1);
    copy_rvec(size,    grid->size);

    return nc_max;
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
            grid->atom_density > 0 ?
            grid->na_c/(grid->atom_density*grid->sx*grid->sy) : 0.0,
            ba[XX], ba[YY], ba[ZZ],
            ba[XX]/grid->sx,
            ba[YY]/grid->sy,
            grid->atom_density > 0 ?
            ba[ZZ]/(grid->na_c/(grid->atom_density*grid->sx*grid->sy)) : 0.0);
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
            grid->atom_density > 0 ?
            grid->na_sc/(grid->atom_density*grid->sx*grid->sy*GPU_NSUBCELL_Z) : 0.0,
            ba[XX], ba[YY], ba[ZZ],
            ba[XX]*GPU_NSUBCELL_X/grid->sx,
            ba[YY]*GPU_NSUBCELL_Y/grid->sy,
            grid->atom_density > 0 ?
            ba[ZZ]/(grid->na_sc/(grid->atom_density*grid->sx*grid->sy*GPU_NSUBCELL_Z)) : 0.0);
}

/* Potentially sorts atoms on LJ coefficients !=0 and ==0.
 * Also sets interaction flags.
 */
void sort_on_lj(int na_c,
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
                /* Only sort when strictly necessary. Ordering particles
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
        sort_on_lj(grid->na_c, a0, a1, atinfo, nbs->a,
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

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ, FALSE, dd_zone,
                   nbs->a+ash, na, x,
                   grid->c0[ZZ],
                   1.0/grid->size[ZZ], ncz*grid->na_sc,
                   sort_work);

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

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ, FALSE, dd_zone,
                   nbs->a+ash, na, x,
                   grid->c0[ZZ],
                   1.0/grid->size[ZZ], ncz*grid->na_sc,
                   sort_work);

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
                grid->bbcz[c*NNBSBB_D  ] = x[nbs->a[ash_z]][ZZ];
                grid->bbcz[c*NNBSBB_D+1] = x[nbs->a[ash_z+na_c-1]][ZZ];
            }

#if GPU_NSUBCELL_Y > 1
            /* Sort the atoms along y */
            sort_atoms(YY, (sub_z & 1), dd_zone,
                       nbs->a+ash_z, na_z, x,
                       grid->c0[YY]+cy*grid->sy,
                       grid->inv_sy, subdiv_z,
                       sort_work);
#endif

            for (sub_y = 0; sub_y < GPU_NSUBCELL_Y; sub_y++)
            {
                ash_y = ash_z + sub_y*subdiv_y;
                na_y  = min(subdiv_y, na-(ash_y-ash));

#if GPU_NSUBCELL_X > 1
                /* Sort the atoms along x */
                sort_atoms(XX, ((cz*GPU_NSUBCELL_Y + sub_y) & 1), dd_zone,
                           nbs->a+ash_y, na_y, x,
                           grid->c0[XX]+cx*grid->sx,
                           grid->inv_sx, subdiv_y,
                           sort_work);
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
    int   cx, cy, cxy, ncz_max, ncz;
    int   nthread, thread;
    int  *cxy_na, cxy_na_i;

    nthread = gmx_omp_nthreads_get(emntPairsearch);

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
        if (nbat->XFormat == nbatX8)
        {
            /* Make the number of cell a multiple of 2 */
            ncz = (ncz + 1) & ~1;
        }
        grid->cxy_ind[i+1] = grid->cxy_ind[i] + ncz;
        /* Clear cxy_na, so we can reuse the array below */
        grid->cxy_na[i] = 0;
    }
    grid->nc = grid->cxy_ind[grid->ncx*grid->ncy] - grid->cxy_ind[0];

    nbat->natoms = (grid->cell0 + grid->nc)*grid->na_sc;

    if (debug)
    {
        fprintf(debug, "ns na_sc %d na_c %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                grid->na_sc, grid->na_c, grid->nc,
                grid->ncx, grid->ncy, grid->nc/((double)(grid->ncx*grid->ncy)),
                ncz_max);
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

    if (dd_zone == 0)
    {
        /* Set the cell indices for the moved particles */
        n0 = grid->nc*grid->na_sc;
        n1 = grid->nc*grid->na_sc+grid->cxy_na[grid->ncx*grid->ncy];
        if (dd_zone == 0)
        {
            for (i = n0; i < n1; i++)
            {
                nbs->cell[nbs->a[i]] = i;
            }
        }
    }

    /* Sort the super-cell columns along z into the sub-cells. */
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

    if (grid->bSimple && nbat->XFormat == nbatX8)
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

    grid->na_c      = nbnxn_kernel_to_ci_size(nb_kernel_type);
    grid->na_cj     = nbnxn_kernel_to_cj_size(nb_kernel_type);
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
        nbs->natoms_nonlocal = max(nbs->natoms_nonlocal, a1);
    }

    /* We always use the home zone (grid[0]) for setting the cell size,
     * since determining densities for non-local zones is difficult.
     */
    nc_max_grid = set_grid_size_xy(nbs, grid,
                                   dd_zone, n-nmoved, corner0, corner1,
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
    int  zone, d;
    rvec c0, c1;

    for (zone = 1; zone < zones->n; zone++)
    {
        for (d = 0; d < DIM; d++)
        {
            c0[d] = zones->size[zone].bb_x0[d];
            c1[d] = zones->size[zone].bb_x1[d];
        }

        nbnxn_put_on_grid(nbs, nbs->ePBC, NULL,
                          zone, c0, c1,
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

    if (grid->bSimple && nbat->XFormat == nbatX8)
    {
        combine_bounding_box_pairs(grid, grid->bb_simple);
    }
}

void nbnxn_get_ncells(nbnxn_search_t nbs, int *ncx, int *ncy)
{
    *ncx = nbs->grid[0].ncx;
    *ncy = nbs->grid[0].ncy;
}

void nbnxn_get_atomorder(nbnxn_search_t nbs, int **a, int *n)
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

/* Determines the cell range along one dimension that
 * the bounding box b0 - b1 sees.
 */
static void get_cell_range(real b0, real b1,
                           int nc, real c0, real s, real invs,
                           real d2, real r2, int *cf, int *cl)
{
    *cf = max((int)((b0 - c0)*invs), 0);

    while (*cf > 0 && d2 + sqr((b0 - c0) - (*cf-1+1)*s) < r2)
    {
        (*cf)--;
    }

    *cl = min((int)((b1 - c0)*invs), nc-1);
    while (*cl < nc-1 && d2 + sqr((*cl+1)*s - (b1 - c0)) < r2)
    {
        (*cl)++;
    }
}

/* Reference code calculating the distance^2 between two bounding boxes */
static float box_dist2(float bx0, float bx1, float by0,
                       float by1, float bz0, float bz1,
                       const nbnxn_bb_t *bb)
{
    float d2;
    float dl, dh, dm, dm0;

    d2 = 0;

    dl  = bx0 - bb->upper[BB_X];
    dh  = bb->lower[BB_X] - bx1;
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    dl  = by0 - bb->upper[BB_Y];
    dh  = bb->lower[BB_Y] - by1;
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    dl  = bz0 - bb->upper[BB_Z];
    dh  = bb->lower[BB_Z] - bz1;
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    return d2;
}

/* Plain C code calculating the distance^2 between two bounding boxes */
static float subc_bb_dist2(int si, const nbnxn_bb_t *bb_i_ci,
                           int csj, const nbnxn_bb_t *bb_j_all)
{
    const nbnxn_bb_t *bb_i, *bb_j;
    float             d2;
    float             dl, dh, dm, dm0;

    bb_i = bb_i_ci  +  si;
    bb_j = bb_j_all + csj;

    d2 = 0;

    dl  = bb_i->lower[BB_X] - bb_j->upper[BB_X];
    dh  = bb_j->lower[BB_X] - bb_i->upper[BB_X];
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    dl  = bb_i->lower[BB_Y] - bb_j->upper[BB_Y];
    dh  = bb_j->lower[BB_Y] - bb_i->upper[BB_Y];
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    dl  = bb_i->lower[BB_Z] - bb_j->upper[BB_Z];
    dh  = bb_j->lower[BB_Z] - bb_i->upper[BB_Z];
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    return d2;
}

#ifdef NBNXN_SEARCH_BB_SIMD4

/* 4-wide SIMD code for bb distance for bb format xyz0 */
static float subc_bb_dist2_simd4(int si, const nbnxn_bb_t *bb_i_ci,
                                 int csj, const nbnxn_bb_t *bb_j_all)
{
    gmx_simd4_float_t bb_i_S0, bb_i_S1;
    gmx_simd4_float_t bb_j_S0, bb_j_S1;
    gmx_simd4_float_t dl_S;
    gmx_simd4_float_t dh_S;
    gmx_simd4_float_t dm_S;
    gmx_simd4_float_t dm0_S;

    bb_i_S0 = gmx_simd4_load_f(&bb_i_ci[si].lower[0]);
    bb_i_S1 = gmx_simd4_load_f(&bb_i_ci[si].upper[0]);
    bb_j_S0 = gmx_simd4_load_f(&bb_j_all[csj].lower[0]);
    bb_j_S1 = gmx_simd4_load_f(&bb_j_all[csj].upper[0]);

    dl_S    = gmx_simd4_sub_f(bb_i_S0, bb_j_S1);
    dh_S    = gmx_simd4_sub_f(bb_j_S0, bb_i_S1);

    dm_S    = gmx_simd4_max_f(dl_S, dh_S);
    dm0_S   = gmx_simd4_max_f(dm_S, gmx_simd4_setzero_f());

    return gmx_simd4_dotproduct3_f(dm0_S, dm0_S);
}

/* Calculate bb bounding distances of bb_i[si,...,si+3] and store them in d2 */
#define SUBC_BB_DIST2_SIMD4_XXXX_INNER(si, bb_i, d2) \
    {                                                \
        int               shi;                                  \
                                                 \
        gmx_simd4_float_t dx_0, dy_0, dz_0;                    \
        gmx_simd4_float_t dx_1, dy_1, dz_1;                    \
                                                 \
        gmx_simd4_float_t mx, my, mz;                          \
        gmx_simd4_float_t m0x, m0y, m0z;                       \
                                                 \
        gmx_simd4_float_t d2x, d2y, d2z;                       \
        gmx_simd4_float_t d2s, d2t;                            \
                                                 \
        shi = si*NNBSBB_D*DIM;                       \
                                                 \
        xi_l = gmx_simd4_load_f(bb_i+shi+0*STRIDE_PBB);   \
        yi_l = gmx_simd4_load_f(bb_i+shi+1*STRIDE_PBB);   \
        zi_l = gmx_simd4_load_f(bb_i+shi+2*STRIDE_PBB);   \
        xi_h = gmx_simd4_load_f(bb_i+shi+3*STRIDE_PBB);   \
        yi_h = gmx_simd4_load_f(bb_i+shi+4*STRIDE_PBB);   \
        zi_h = gmx_simd4_load_f(bb_i+shi+5*STRIDE_PBB);   \
                                                 \
        dx_0 = gmx_simd4_sub_f(xi_l, xj_h);                 \
        dy_0 = gmx_simd4_sub_f(yi_l, yj_h);                 \
        dz_0 = gmx_simd4_sub_f(zi_l, zj_h);                 \
                                                 \
        dx_1 = gmx_simd4_sub_f(xj_l, xi_h);                 \
        dy_1 = gmx_simd4_sub_f(yj_l, yi_h);                 \
        dz_1 = gmx_simd4_sub_f(zj_l, zi_h);                 \
                                                 \
        mx   = gmx_simd4_max_f(dx_0, dx_1);                 \
        my   = gmx_simd4_max_f(dy_0, dy_1);                 \
        mz   = gmx_simd4_max_f(dz_0, dz_1);                 \
                                                 \
        m0x  = gmx_simd4_max_f(mx, zero);                   \
        m0y  = gmx_simd4_max_f(my, zero);                   \
        m0z  = gmx_simd4_max_f(mz, zero);                   \
                                                 \
        d2x  = gmx_simd4_mul_f(m0x, m0x);                   \
        d2y  = gmx_simd4_mul_f(m0y, m0y);                   \
        d2z  = gmx_simd4_mul_f(m0z, m0z);                   \
                                                 \
        d2s  = gmx_simd4_add_f(d2x, d2y);                   \
        d2t  = gmx_simd4_add_f(d2s, d2z);                   \
                                                 \
        gmx_simd4_store_f(d2+si, d2t);                      \
    }

/* 4-wide SIMD code for nsi bb distances for bb format xxxxyyyyzzzz */
static void subc_bb_dist2_simd4_xxxx(const float *bb_j,
                                     int nsi, const float *bb_i,
                                     float *d2)
{
    gmx_simd4_float_t xj_l, yj_l, zj_l;
    gmx_simd4_float_t xj_h, yj_h, zj_h;
    gmx_simd4_float_t xi_l, yi_l, zi_l;
    gmx_simd4_float_t xi_h, yi_h, zi_h;

    gmx_simd4_float_t zero;

    zero = gmx_simd4_setzero_f();

    xj_l = gmx_simd4_set1_f(bb_j[0*STRIDE_PBB]);
    yj_l = gmx_simd4_set1_f(bb_j[1*STRIDE_PBB]);
    zj_l = gmx_simd4_set1_f(bb_j[2*STRIDE_PBB]);
    xj_h = gmx_simd4_set1_f(bb_j[3*STRIDE_PBB]);
    yj_h = gmx_simd4_set1_f(bb_j[4*STRIDE_PBB]);
    zj_h = gmx_simd4_set1_f(bb_j[5*STRIDE_PBB]);

    /* Here we "loop" over si (0,STRIDE_PBB) from 0 to nsi with step STRIDE_PBB.
     * But as we know the number of iterations is 1 or 2, we unroll manually.
     */
    SUBC_BB_DIST2_SIMD4_XXXX_INNER(0, bb_i, d2);
    if (STRIDE_PBB < nsi)
    {
        SUBC_BB_DIST2_SIMD4_XXXX_INNER(STRIDE_PBB, bb_i, d2);
    }
}

#endif /* NBNXN_SEARCH_BB_SIMD4 */

/* Plain C function which determines if any atom pair between two cells
 * is within distance sqrt(rl2).
 */
static gmx_bool subc_in_range_x(int na_c,
                                int si, const real *x_i,
                                int csj, int stride, const real *x_j,
                                real rl2)
{
    int  i, j, i0, j0;
    real d2;

    for (i = 0; i < na_c; i++)
    {
        i0 = (si*na_c + i)*DIM;
        for (j = 0; j < na_c; j++)
        {
            j0 = (csj*na_c + j)*stride;

            d2 = sqr(x_i[i0  ] - x_j[j0  ]) +
                sqr(x_i[i0+1] - x_j[j0+1]) +
                sqr(x_i[i0+2] - x_j[j0+2]);

            if (d2 < rl2)
            {
                return TRUE;
            }
        }
    }

    return FALSE;
}

#ifdef NBNXN_SEARCH_SIMD4_FLOAT_X_BB

/* 4-wide SIMD function which determines if any atom pair between two cells,
 * both with 8 atoms, is within distance sqrt(rl2).
 * Using 8-wide AVX is not faster on Intel Sandy Bridge.
 */
static gmx_bool subc_in_range_simd4(int na_c,
                                    int si, const real *x_i,
                                    int csj, int stride, const real *x_j,
                                    real rl2)
{
    gmx_simd4_real_t ix_S0, iy_S0, iz_S0;
    gmx_simd4_real_t ix_S1, iy_S1, iz_S1;

    gmx_simd4_real_t rc2_S;

    int              dim_stride;
    int              j0, j1;

    rc2_S   = gmx_simd4_set1_r(rl2);

    dim_stride = NBNXN_GPU_CLUSTER_SIZE/STRIDE_PBB*DIM;
    ix_S0      = gmx_simd4_load_r(x_i+(si*dim_stride+0)*STRIDE_PBB);
    iy_S0      = gmx_simd4_load_r(x_i+(si*dim_stride+1)*STRIDE_PBB);
    iz_S0      = gmx_simd4_load_r(x_i+(si*dim_stride+2)*STRIDE_PBB);
    ix_S1      = gmx_simd4_load_r(x_i+(si*dim_stride+3)*STRIDE_PBB);
    iy_S1      = gmx_simd4_load_r(x_i+(si*dim_stride+4)*STRIDE_PBB);
    iz_S1      = gmx_simd4_load_r(x_i+(si*dim_stride+5)*STRIDE_PBB);

    /* We loop from the outer to the inner particles to maximize
     * the chance that we find a pair in range quickly and return.
     */
    j0 = csj*na_c;
    j1 = j0 + na_c - 1;
    while (j0 < j1)
    {
        gmx_simd4_real_t jx0_S, jy0_S, jz0_S;
        gmx_simd4_real_t jx1_S, jy1_S, jz1_S;

        gmx_simd4_real_t dx_S0, dy_S0, dz_S0;
        gmx_simd4_real_t dx_S1, dy_S1, dz_S1;
        gmx_simd4_real_t dx_S2, dy_S2, dz_S2;
        gmx_simd4_real_t dx_S3, dy_S3, dz_S3;

        gmx_simd4_real_t rsq_S0;
        gmx_simd4_real_t rsq_S1;
        gmx_simd4_real_t rsq_S2;
        gmx_simd4_real_t rsq_S3;

        gmx_simd4_bool_t wco_S0;
        gmx_simd4_bool_t wco_S1;
        gmx_simd4_bool_t wco_S2;
        gmx_simd4_bool_t wco_S3;
        gmx_simd4_bool_t wco_any_S01, wco_any_S23, wco_any_S;

        jx0_S = gmx_simd4_set1_r(x_j[j0*stride+0]);
        jy0_S = gmx_simd4_set1_r(x_j[j0*stride+1]);
        jz0_S = gmx_simd4_set1_r(x_j[j0*stride+2]);

        jx1_S = gmx_simd4_set1_r(x_j[j1*stride+0]);
        jy1_S = gmx_simd4_set1_r(x_j[j1*stride+1]);
        jz1_S = gmx_simd4_set1_r(x_j[j1*stride+2]);

        /* Calculate distance */
        dx_S0            = gmx_simd4_sub_r(ix_S0, jx0_S);
        dy_S0            = gmx_simd4_sub_r(iy_S0, jy0_S);
        dz_S0            = gmx_simd4_sub_r(iz_S0, jz0_S);
        dx_S1            = gmx_simd4_sub_r(ix_S1, jx0_S);
        dy_S1            = gmx_simd4_sub_r(iy_S1, jy0_S);
        dz_S1            = gmx_simd4_sub_r(iz_S1, jz0_S);
        dx_S2            = gmx_simd4_sub_r(ix_S0, jx1_S);
        dy_S2            = gmx_simd4_sub_r(iy_S0, jy1_S);
        dz_S2            = gmx_simd4_sub_r(iz_S0, jz1_S);
        dx_S3            = gmx_simd4_sub_r(ix_S1, jx1_S);
        dy_S3            = gmx_simd4_sub_r(iy_S1, jy1_S);
        dz_S3            = gmx_simd4_sub_r(iz_S1, jz1_S);

        /* rsq = dx*dx+dy*dy+dz*dz */
        rsq_S0           = gmx_simd4_calc_rsq_r(dx_S0, dy_S0, dz_S0);
        rsq_S1           = gmx_simd4_calc_rsq_r(dx_S1, dy_S1, dz_S1);
        rsq_S2           = gmx_simd4_calc_rsq_r(dx_S2, dy_S2, dz_S2);
        rsq_S3           = gmx_simd4_calc_rsq_r(dx_S3, dy_S3, dz_S3);

        wco_S0           = gmx_simd4_cmplt_r(rsq_S0, rc2_S);
        wco_S1           = gmx_simd4_cmplt_r(rsq_S1, rc2_S);
        wco_S2           = gmx_simd4_cmplt_r(rsq_S2, rc2_S);
        wco_S3           = gmx_simd4_cmplt_r(rsq_S3, rc2_S);

        wco_any_S01      = gmx_simd4_or_b(wco_S0, wco_S1);
        wco_any_S23      = gmx_simd4_or_b(wco_S2, wco_S3);
        wco_any_S        = gmx_simd4_or_b(wco_any_S01, wco_any_S23);

        if (gmx_simd4_anytrue_b(wco_any_S))
        {
            return TRUE;
        }

        j0++;
        j1--;
    }
    return FALSE;

}
#endif


/* Returns the j sub-cell for index cj_ind */
static int nbl_cj(const nbnxn_pairlist_t *nbl, int cj_ind)
{
    return nbl->cj4[cj_ind >> NBNXN_GPU_JGROUP_SIZE_2LOG].cj[cj_ind & (NBNXN_GPU_JGROUP_SIZE - 1)];
}

/* Returns the i-interaction mask of the j sub-cell for index cj_ind */
static unsigned int nbl_imask0(const nbnxn_pairlist_t *nbl, int cj_ind)
{
    return nbl->cj4[cj_ind >> NBNXN_GPU_JGROUP_SIZE_2LOG].imei[0].imask;
}

/* Ensures there is enough space for extra extra exclusion masks */
static void check_excl_space(nbnxn_pairlist_t *nbl, int extra)
{
    if (nbl->nexcl+extra > nbl->excl_nalloc)
    {
        nbl->excl_nalloc = over_alloc_small(nbl->nexcl+extra);
        nbnxn_realloc_void((void **)&nbl->excl,
                           nbl->nexcl*sizeof(*nbl->excl),
                           nbl->excl_nalloc*sizeof(*nbl->excl),
                           nbl->alloc, nbl->free);
    }
}

/* Ensures there is enough space for ncell extra j-cells in the list */
static void check_subcell_list_space_simple(nbnxn_pairlist_t *nbl,
                                            int               ncell)
{
    int cj_max;

    cj_max = nbl->ncj + ncell;

    if (cj_max > nbl->cj_nalloc)
    {
        nbl->cj_nalloc = over_alloc_small(cj_max);
        nbnxn_realloc_void((void **)&nbl->cj,
                           nbl->ncj*sizeof(*nbl->cj),
                           nbl->cj_nalloc*sizeof(*nbl->cj),
                           nbl->alloc, nbl->free);
    }
}

/* Ensures there is enough space for ncell extra j-subcells in the list */
static void check_subcell_list_space_supersub(nbnxn_pairlist_t *nbl,
                                              int               nsupercell)
{
    int ncj4_max, j4, j, w, t;

#define NWARP       2
#define WARP_SIZE  32

    /* We can have maximally nsupercell*GPU_NSUBCELL sj lists */
    /* We can store 4 j-subcell - i-supercell pairs in one struct.
     * since we round down, we need one extra entry.
     */
    ncj4_max = ((nbl->work->cj_ind + nsupercell*GPU_NSUBCELL + NBNXN_GPU_JGROUP_SIZE - 1) >> NBNXN_GPU_JGROUP_SIZE_2LOG);

    if (ncj4_max > nbl->cj4_nalloc)
    {
        nbl->cj4_nalloc = over_alloc_small(ncj4_max);
        nbnxn_realloc_void((void **)&nbl->cj4,
                           nbl->work->cj4_init*sizeof(*nbl->cj4),
                           nbl->cj4_nalloc*sizeof(*nbl->cj4),
                           nbl->alloc, nbl->free);
    }

    if (ncj4_max > nbl->work->cj4_init)
    {
        for (j4 = nbl->work->cj4_init; j4 < ncj4_max; j4++)
        {
            /* No i-subcells and no excl's in the list initially */
            for (w = 0; w < NWARP; w++)
            {
                nbl->cj4[j4].imei[w].imask    = 0U;
                nbl->cj4[j4].imei[w].excl_ind = 0;

            }
        }
        nbl->work->cj4_init = ncj4_max;
    }
}

/* Set all excl masks for one GPU warp no exclusions */
static void set_no_excls(nbnxn_excl_t *excl)
{
    int t;

    for (t = 0; t < WARP_SIZE; t++)
    {
        /* Turn all interaction bits on */
        excl->pair[t] = NBNXN_INTERACTION_MASK_ALL;
    }
}

/* Initializes a single nbnxn_pairlist_t data structure */
static void nbnxn_init_pairlist(nbnxn_pairlist_t *nbl,
                                gmx_bool          bSimple,
                                nbnxn_alloc_t    *alloc,
                                nbnxn_free_t     *free)
{
    if (alloc == NULL)
    {
        nbl->alloc = nbnxn_alloc_aligned;
    }
    else
    {
        nbl->alloc = alloc;
    }
    if (free == NULL)
    {
        nbl->free = nbnxn_free_aligned;
    }
    else
    {
        nbl->free = free;
    }

    nbl->bSimple     = bSimple;
    nbl->na_sc       = 0;
    nbl->na_ci       = 0;
    nbl->na_cj       = 0;
    nbl->nci         = 0;
    nbl->ci          = NULL;
    nbl->ci_nalloc   = 0;
    nbl->ncj         = 0;
    nbl->cj          = NULL;
    nbl->cj_nalloc   = 0;
    nbl->ncj4        = 0;
    /* We need one element extra in sj, so alloc initially with 1 */
    nbl->cj4_nalloc  = 0;
    nbl->cj4         = NULL;
    nbl->nci_tot     = 0;

    if (!nbl->bSimple)
    {
        nbl->excl        = NULL;
        nbl->excl_nalloc = 0;
        nbl->nexcl       = 0;
        check_excl_space(nbl, 1);
        nbl->nexcl       = 1;
        set_no_excls(&nbl->excl[0]);
    }

    snew(nbl->work, 1);
    if (nbl->bSimple)
    {
        snew_aligned(nbl->work->bb_ci, 1, NBNXN_SEARCH_BB_MEM_ALIGN);
    }
    else
    {
#ifdef NBNXN_BBXXXX
        snew_aligned(nbl->work->pbb_ci, GPU_NSUBCELL/STRIDE_PBB*NNBSBB_XXXX, NBNXN_SEARCH_BB_MEM_ALIGN);
#else
        snew_aligned(nbl->work->bb_ci, GPU_NSUBCELL, NBNXN_SEARCH_BB_MEM_ALIGN);
#endif
    }
    snew_aligned(nbl->work->x_ci, NBNXN_NA_SC_MAX*DIM, NBNXN_SEARCH_BB_MEM_ALIGN);
#ifdef GMX_NBNXN_SIMD
    snew_aligned(nbl->work->x_ci_simd_4xn, 1, NBNXN_MEM_ALIGN);
    snew_aligned(nbl->work->x_ci_simd_2xnn, 1, NBNXN_MEM_ALIGN);
#endif
    snew_aligned(nbl->work->d2, GPU_NSUBCELL, NBNXN_SEARCH_BB_MEM_ALIGN);

    nbl->work->sort            = NULL;
    nbl->work->sort_nalloc     = 0;
    nbl->work->sci_sort        = NULL;
    nbl->work->sci_sort_nalloc = 0;
}

void nbnxn_init_pairlist_set(nbnxn_pairlist_set_t *nbl_list,
                             gmx_bool bSimple, gmx_bool bCombined,
                             nbnxn_alloc_t *alloc,
                             nbnxn_free_t  *free)
{
    int i;

    nbl_list->bSimple   = bSimple;
    nbl_list->bCombined = bCombined;

    nbl_list->nnbl = gmx_omp_nthreads_get(emntNonbonded);

    if (!nbl_list->bCombined &&
        nbl_list->nnbl > NBNXN_BUFFERFLAG_MAX_THREADS)
    {
        gmx_fatal(FARGS, "%d OpenMP threads were requested. Since the non-bonded force buffer reduction is prohibitively slow with more than %d threads, we do not allow this. Use %d or less OpenMP threads.",
                  nbl_list->nnbl, NBNXN_BUFFERFLAG_MAX_THREADS, NBNXN_BUFFERFLAG_MAX_THREADS);
    }

    snew(nbl_list->nbl, nbl_list->nnbl);
    snew(nbl_list->nbl_fep, nbl_list->nnbl);
    /* Execute in order to avoid memory interleaving between threads */
#pragma omp parallel for num_threads(nbl_list->nnbl) schedule(static)
    for (i = 0; i < nbl_list->nnbl; i++)
    {
        /* Allocate the nblist data structure locally on each thread
         * to optimize memory access for NUMA architectures.
         */
        snew(nbl_list->nbl[i], 1);

        /* Only list 0 is used on the GPU, use normal allocation for i>0 */
        if (i == 0)
        {
            nbnxn_init_pairlist(nbl_list->nbl[i], nbl_list->bSimple, alloc, free);
        }
        else
        {
            nbnxn_init_pairlist(nbl_list->nbl[i], nbl_list->bSimple, NULL, NULL);
        }

        snew(nbl_list->nbl_fep[i], 1);
        nbnxn_init_pairlist_fep(nbl_list->nbl_fep[i]);
    }
}

/* Print statistics of a pair list, used for debug output */
static void print_nblist_statistics_simple(FILE *fp, const nbnxn_pairlist_t *nbl,
                                           const nbnxn_search_t nbs, real rl)
{
    const nbnxn_grid_t *grid;
    int                 cs[SHIFTS];
    int                 s, i, j;
    int                 npexcl;

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp, "nbl nci %d ncj %d\n",
            nbl->nci, nbl->ncj);
    fprintf(fp, "nbl na_sc %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->na_sc, rl, nbl->ncj, nbl->ncj/(double)grid->nc,
            nbl->ncj/(double)grid->nc*grid->na_sc,
            nbl->ncj/(double)grid->nc*grid->na_sc/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nc*grid->na_sc/(grid->size[XX]*grid->size[YY]*grid->size[ZZ])));

    fprintf(fp, "nbl average j cell list length %.1f\n",
            0.25*nbl->ncj/(double)max(nbl->nci, 1));

    for (s = 0; s < SHIFTS; s++)
    {
        cs[s] = 0;
    }
    npexcl = 0;
    for (i = 0; i < nbl->nci; i++)
    {
        cs[nbl->ci[i].shift & NBNXN_CI_SHIFT] +=
            nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start;

        j = nbl->ci[i].cj_ind_start;
        while (j < nbl->ci[i].cj_ind_end &&
               nbl->cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            npexcl++;
            j++;
        }
    }
    fprintf(fp, "nbl cell pairs, total: %d excl: %d %.1f%%\n",
            nbl->ncj, npexcl, 100*npexcl/(double)max(nbl->ncj, 1));
    for (s = 0; s < SHIFTS; s++)
    {
        if (cs[s] > 0)
        {
            fprintf(fp, "nbl shift %2d ncj %3d\n", s, cs[s]);
        }
    }
}

/* Print statistics of a pair lists, used for debug output */
static void print_nblist_statistics_supersub(FILE *fp, const nbnxn_pairlist_t *nbl,
                                             const nbnxn_search_t nbs, real rl)
{
    const nbnxn_grid_t *grid;
    int                 i, j4, j, si, b;
    int                 c[GPU_NSUBCELL+1];
    double              sum_nsp, sum_nsp2;
    int                 nsp_max;

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp, "nbl nsci %d ncj4 %d nsi %d excl4 %d\n",
            nbl->nsci, nbl->ncj4, nbl->nci_tot, nbl->nexcl);
    fprintf(fp, "nbl na_c %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->na_ci, rl, nbl->nci_tot, nbl->nci_tot/(double)grid->nsubc_tot,
            nbl->nci_tot/(double)grid->nsubc_tot*grid->na_c,
            nbl->nci_tot/(double)grid->nsubc_tot*grid->na_c/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nsubc_tot*grid->na_c/(grid->size[XX]*grid->size[YY]*grid->size[ZZ])));

    sum_nsp  = 0;
    sum_nsp2 = 0;
    nsp_max  = 0;
    for (si = 0; si <= GPU_NSUBCELL; si++)
    {
        c[si] = 0;
    }
    for (i = 0; i < nbl->nsci; i++)
    {
        int nsp;

        nsp = 0;
        for (j4 = nbl->sci[i].cj4_ind_start; j4 < nbl->sci[i].cj4_ind_end; j4++)
        {
            for (j = 0; j < NBNXN_GPU_JGROUP_SIZE; j++)
            {
                b = 0;
                for (si = 0; si < GPU_NSUBCELL; si++)
                {
                    if (nbl->cj4[j4].imei[0].imask & (1U << (j*GPU_NSUBCELL + si)))
                    {
                        b++;
                    }
                }
                nsp += b;
                c[b]++;
            }
        }
        sum_nsp  += nsp;
        sum_nsp2 += nsp*nsp;
        nsp_max   = max(nsp_max, nsp);
    }
    if (nbl->nsci > 0)
    {
        sum_nsp  /= nbl->nsci;
        sum_nsp2 /= nbl->nsci;
    }
    fprintf(fp, "nbl #cluster-pairs: av %.1f stddev %.1f max %d\n",
            sum_nsp, sqrt(sum_nsp2 - sum_nsp*sum_nsp), nsp_max);

    if (nbl->ncj4 > 0)
    {
        for (b = 0; b <= GPU_NSUBCELL; b++)
        {
            fprintf(fp, "nbl j-list #i-subcell %d %7d %4.1f\n",
                    b, c[b],
                    100.0*c[b]/(double)(nbl->ncj4*NBNXN_GPU_JGROUP_SIZE));
        }
    }
}

/* Returns a pointer to the exclusion mask for cj4-unit cj4, warp warp */
static void low_get_nbl_exclusions(nbnxn_pairlist_t *nbl, int cj4,
                                   int warp, nbnxn_excl_t **excl)
{
    if (nbl->cj4[cj4].imei[warp].excl_ind == 0)
    {
        /* No exclusions set, make a new list entry */
        nbl->cj4[cj4].imei[warp].excl_ind = nbl->nexcl;
        nbl->nexcl++;
        *excl = &nbl->excl[nbl->cj4[cj4].imei[warp].excl_ind];
        set_no_excls(*excl);
    }
    else
    {
        /* We already have some exclusions, new ones can be added to the list */
        *excl = &nbl->excl[nbl->cj4[cj4].imei[warp].excl_ind];
    }
}

/* Returns a pointer to the exclusion mask for cj4-unit cj4, warp warp,
 * generates a new element and allocates extra memory, if necessary.
 */
static void get_nbl_exclusions_1(nbnxn_pairlist_t *nbl, int cj4,
                                 int warp, nbnxn_excl_t **excl)
{
    if (nbl->cj4[cj4].imei[warp].excl_ind == 0)
    {
        /* We need to make a new list entry, check if we have space */
        check_excl_space(nbl, 1);
    }
    low_get_nbl_exclusions(nbl, cj4, warp, excl);
}

/* Returns pointers to the exclusion mask for cj4-unit cj4 for both warps,
 * generates a new element and allocates extra memory, if necessary.
 */
static void get_nbl_exclusions_2(nbnxn_pairlist_t *nbl, int cj4,
                                 nbnxn_excl_t **excl_w0,
                                 nbnxn_excl_t **excl_w1)
{
    /* Check for space we might need */
    check_excl_space(nbl, 2);

    low_get_nbl_exclusions(nbl, cj4, 0, excl_w0);
    low_get_nbl_exclusions(nbl, cj4, 1, excl_w1);
}

/* Sets the self exclusions i=j and pair exclusions i>j */
static void set_self_and_newton_excls_supersub(nbnxn_pairlist_t *nbl,
                                               int cj4_ind, int sj_offset,
                                               int si)
{
    nbnxn_excl_t *excl[2];
    int           ei, ej, w;

    /* Here we only set the set self and double pair exclusions */

    get_nbl_exclusions_2(nbl, cj4_ind, &excl[0], &excl[1]);

    /* Only minor < major bits set */
    for (ej = 0; ej < nbl->na_ci; ej++)
    {
        w = (ej>>2);
        for (ei = ej; ei < nbl->na_ci; ei++)
        {
            excl[w]->pair[(ej & (NBNXN_GPU_JGROUP_SIZE-1))*nbl->na_ci + ei] &=
                ~(1U << (sj_offset*GPU_NSUBCELL + si));
        }
    }
}

/* Returns a diagonal or off-diagonal interaction mask for plain C lists */
static unsigned int get_imask(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=2 */
static unsigned int get_imask_simd_j2(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci*2 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_0 :
            (rdiag && ci*2+1 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_1 :
             NBNXN_INTERACTION_MASK_ALL));
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=4 */
static unsigned int get_imask_simd_j4(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=8 */
static unsigned int get_imask_simd_j8(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj*2 ? NBNXN_INTERACTION_MASK_DIAG_J8_0 :
            (rdiag && ci == cj*2+1 ? NBNXN_INTERACTION_MASK_DIAG_J8_1 :
             NBNXN_INTERACTION_MASK_ALL));
}

#ifdef GMX_NBNXN_SIMD
#if GMX_SIMD_REAL_WIDTH == 2
#define get_imask_simd_4xn  get_imask_simd_j2
#endif
#if GMX_SIMD_REAL_WIDTH == 4
#define get_imask_simd_4xn  get_imask_simd_j4
#endif
#if GMX_SIMD_REAL_WIDTH == 8
#define get_imask_simd_4xn  get_imask_simd_j8
#define get_imask_simd_2xnn get_imask_simd_j4
#endif
#if GMX_SIMD_REAL_WIDTH == 16
#define get_imask_simd_2xnn get_imask_simd_j8
#endif
#endif

/* Plain C code for making a pair list of cell ci vs cell cjf-cjl.
 * Checks bounding box distances and possibly atom pair distances.
 */
static void make_cluster_list_simple(const nbnxn_grid_t *gridj,
                                     nbnxn_pairlist_t *nbl,
                                     int ci, int cjf, int cjl,
                                     gmx_bool remove_sub_diag,
                                     const real *x_j,
                                     real rl2, float rbb2,
                                     int *ndistc)
{
    const nbnxn_list_work_t *work;

    const nbnxn_bb_t        *bb_ci;
    const real              *x_ci;

    gmx_bool                 InRange;
    real                     d2;
    int                      cjf_gl, cjl_gl, cj;

    work = nbl->work;

    bb_ci = nbl->work->bb_ci;
    x_ci  = nbl->work->x_ci;

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
        d2       = subc_bb_dist2(0, bb_ci, cjf, gridj->bb);
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
            int i, j;

            cjf_gl = gridj->cell0 + cjf;
            for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                {
                    InRange = InRange ||
                        (sqr(x_ci[i*STRIDE_XYZ+XX] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+XX]) +
                         sqr(x_ci[i*STRIDE_XYZ+YY] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+YY]) +
                         sqr(x_ci[i*STRIDE_XYZ+ZZ] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+ZZ]) < rl2);
                }
            }
            *ndistc += NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
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
        d2       = subc_bb_dist2(0, bb_ci, cjl, gridj->bb);
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
            int i, j;

            cjl_gl = gridj->cell0 + cjl;
            for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                {
                    InRange = InRange ||
                        (sqr(x_ci[i*STRIDE_XYZ+XX] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+XX]) +
                         sqr(x_ci[i*STRIDE_XYZ+YY] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+YY]) +
                         sqr(x_ci[i*STRIDE_XYZ+ZZ] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+ZZ]) < rl2);
                }
            }
            *ndistc += NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
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
            nbl->cj[nbl->ncj].cj   = gridj->cell0 + cj;
            nbl->cj[nbl->ncj].excl = get_imask(remove_sub_diag, ci, cj);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/mdlib/nbnxn_search_simd_4xn.h"
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/mdlib/nbnxn_search_simd_2xnn.h"
#endif

/* Plain C or SIMD4 code for making a pair list of super-cell sci vs scj.
 * Checks bounding box distances and possibly atom pair distances.
 */
static void make_cluster_list_supersub(const nbnxn_grid_t *gridi,
                                       const nbnxn_grid_t *gridj,
                                       nbnxn_pairlist_t *nbl,
                                       int sci, int scj,
                                       gmx_bool sci_equals_scj,
                                       int stride, const real *x,
                                       real rl2, float rbb2,
                                       int *ndistc)
{
    int               na_c;
    int               npair;
    int               cjo, ci1, ci, cj, cj_gl;
    int               cj4_ind, cj_offset;
    unsigned int      imask;
    nbnxn_cj4_t      *cj4;
#ifdef NBNXN_BBXXXX
    const float      *pbb_ci;
#else
    const nbnxn_bb_t *bb_ci;
#endif
    const real       *x_ci;
    float            *d2l, d2;
    int               w;
#define PRUNE_LIST_CPU_ONE
#ifdef PRUNE_LIST_CPU_ONE
    int  ci_last = -1;
#endif

    d2l = nbl->work->d2;

#ifdef NBNXN_BBXXXX
    pbb_ci = nbl->work->pbb_ci;
#else
    bb_ci  = nbl->work->bb_ci;
#endif
    x_ci   = nbl->work->x_ci;

    na_c = gridj->na_c;

    for (cjo = 0; cjo < gridj->nsubc[scj]; cjo++)
    {
        cj4_ind   = (nbl->work->cj_ind >> NBNXN_GPU_JGROUP_SIZE_2LOG);
        cj_offset = nbl->work->cj_ind - cj4_ind*NBNXN_GPU_JGROUP_SIZE;
        cj4       = &nbl->cj4[cj4_ind];

        cj = scj*GPU_NSUBCELL + cjo;

        cj_gl = gridj->cell0*GPU_NSUBCELL + cj;

        /* Initialize this j-subcell i-subcell list */
        cj4->cj[cj_offset] = cj_gl;
        imask              = 0;

        if (sci_equals_scj)
        {
            ci1 = cjo + 1;
        }
        else
        {
            ci1 = gridi->nsubc[sci];
        }

#ifdef NBNXN_BBXXXX
        /* Determine all ci1 bb distances in one call with SIMD4 */
        subc_bb_dist2_simd4_xxxx(gridj->pbb+(cj>>STRIDE_PBB_2LOG)*NNBSBB_XXXX+(cj & (STRIDE_PBB-1)),
                                 ci1, pbb_ci, d2l);
        *ndistc += na_c*2;
#endif

        npair = 0;
        /* We use a fixed upper-bound instead of ci1 to help optimization */
        for (ci = 0; ci < GPU_NSUBCELL; ci++)
        {
            if (ci == ci1)
            {
                break;
            }

#ifndef NBNXN_BBXXXX
            /* Determine the bb distance between ci and cj */
            d2l[ci]  = subc_bb_dist2(ci, bb_ci, cj, gridj->bb);
            *ndistc += 2;
#endif
            d2 = d2l[ci];

#ifdef PRUNE_LIST_CPU_ALL
            /* Check if the distance is within the distance where
             * we use only the bounding box distance rbb,
             * or within the cut-off and there is at least one atom pair
             * within the cut-off. This check is very costly.
             */
            *ndistc += na_c*na_c;
            if (d2 < rbb2 ||
                (d2 < rl2 &&
#ifdef NBNXN_PBB_SIMD4
                 subc_in_range_simd4
#else
                 subc_in_range_x
#endif
                     (na_c, ci, x_ci, cj_gl, stride, x, rl2)))
#else
            /* Check if the distance between the two bounding boxes
             * in within the pair-list cut-off.
             */
            if (d2 < rl2)
#endif
            {
                /* Flag this i-subcell to be taken into account */
                imask |= (1U << (cj_offset*GPU_NSUBCELL+ci));

#ifdef PRUNE_LIST_CPU_ONE
                ci_last = ci;
#endif

                npair++;
            }
        }

#ifdef PRUNE_LIST_CPU_ONE
        /* If we only found 1 pair, check if any atoms are actually
         * within the cut-off, so we could get rid of it.
         */
        if (npair == 1 && d2l[ci_last] >= rbb2)
        {
            /* Avoid using function pointers here, as it's slower */
            if (
#ifdef NBNXN_PBB_SIMD4
                !subc_in_range_simd4
#else
                !subc_in_range_x
#endif
                    (na_c, ci_last, x_ci, cj_gl, stride, x, rl2))
            {
                imask &= ~(1U << (cj_offset*GPU_NSUBCELL+ci_last));
                npair--;
            }
        }
#endif

        if (npair > 0)
        {
            /* We have a useful sj entry, close it now */

            /* Set the exclucions for the ci== sj entry.
             * Here we don't bother to check if this entry is actually flagged,
             * as it will nearly always be in the list.
             */
            if (sci_equals_scj)
            {
                set_self_and_newton_excls_supersub(nbl, cj4_ind, cj_offset, cjo);
            }

            /* Copy the cluster interaction mask to the list */
            for (w = 0; w < NWARP; w++)
            {
                cj4->imei[w].imask |= imask;
            }

            nbl->work->cj_ind++;

            /* Keep the count */
            nbl->nci_tot += npair;

            /* Increase the closing index in i super-cell list */
            nbl->sci[nbl->nsci].cj4_ind_end =
                ((nbl->work->cj_ind+NBNXN_GPU_JGROUP_SIZE-1) >> NBNXN_GPU_JGROUP_SIZE_2LOG);
        }
    }
}

/* Set all atom-pair exclusions from the topology stored in excl
 * as masks in the pair-list for simple list i-entry nbl_ci
 */
static void set_ci_top_excls(const nbnxn_search_t nbs,
                             nbnxn_pairlist_t    *nbl,
                             gmx_bool             diagRemoved,
                             int                  na_ci_2log,
                             int                  na_cj_2log,
                             const nbnxn_ci_t    *nbl_ci,
                             const t_blocka      *excl)
{
    const int    *cell;
    int           ci;
    int           cj_ind_first, cj_ind_last;
    int           cj_first, cj_last;
    int           ndirect;
    int           i, ai, aj, si, eind, ge, se;
    int           found, cj_ind_0, cj_ind_1, cj_ind_m;
    int           cj_m;
    gmx_bool      Found_si;
    int           si_ind;
    nbnxn_excl_t *nbl_excl;
    int           inner_i, inner_e;

    cell = nbs->cell;

    if (nbl_ci->cj_ind_end == nbl_ci->cj_ind_start)
    {
        /* Empty list */
        return;
    }

    ci = nbl_ci->ci;

    cj_ind_first = nbl_ci->cj_ind_start;
    cj_ind_last  = nbl->ncj - 1;

    cj_first = nbl->cj[cj_ind_first].cj;
    cj_last  = nbl->cj[cj_ind_last].cj;

    /* Determine how many contiguous j-cells we have starting
     * from the first i-cell. This number can be used to directly
     * calculate j-cell indices for excluded atoms.
     */
    ndirect = 0;
    if (na_ci_2log == na_cj_2log)
    {
        while (cj_ind_first + ndirect <= cj_ind_last &&
               nbl->cj[cj_ind_first+ndirect].cj == ci + ndirect)
        {
            ndirect++;
        }
    }
#ifdef NBNXN_SEARCH_BB_SIMD4
    else
    {
        while (cj_ind_first + ndirect <= cj_ind_last &&
               nbl->cj[cj_ind_first+ndirect].cj == ci_to_cj(na_cj_2log, ci) + ndirect)
        {
            ndirect++;
        }
    }
#endif

    /* Loop over the atoms in the i super-cell */
    for (i = 0; i < nbl->na_sc; i++)
    {
        ai = nbs->a[ci*nbl->na_sc+i];
        if (ai >= 0)
        {
            si  = (i>>na_ci_2log);

            /* Loop over the topology-based exclusions for this i-atom */
            for (eind = excl->index[ai]; eind < excl->index[ai+1]; eind++)
            {
                aj = excl->a[eind];

                if (aj == ai)
                {
                    /* The self exclusion are already set, save some time */
                    continue;
                }

                ge = cell[aj];

                /* Without shifts we only calculate interactions j>i
                 * for one-way pair-lists.
                 */
                if (diagRemoved && ge <= ci*nbl->na_sc + i)
                {
                    continue;
                }

                se = (ge >> na_cj_2log);

                /* Could the cluster se be in our list? */
                if (se >= cj_first && se <= cj_last)
                {
                    if (se < cj_first + ndirect)
                    {
                        /* We can calculate cj_ind directly from se */
                        found = cj_ind_first + se - cj_first;
                    }
                    else
                    {
                        /* Search for se using bisection */
                        found    = -1;
                        cj_ind_0 = cj_ind_first + ndirect;
                        cj_ind_1 = cj_ind_last + 1;
                        while (found == -1 && cj_ind_0 < cj_ind_1)
                        {
                            cj_ind_m = (cj_ind_0 + cj_ind_1)>>1;

                            cj_m = nbl->cj[cj_ind_m].cj;

                            if (se == cj_m)
                            {
                                found = cj_ind_m;
                            }
                            else if (se < cj_m)
                            {
                                cj_ind_1 = cj_ind_m;
                            }
                            else
                            {
                                cj_ind_0 = cj_ind_m + 1;
                            }
                        }
                    }

                    if (found >= 0)
                    {
                        inner_i = i  - (si << na_ci_2log);
                        inner_e = ge - (se << na_cj_2log);

                        nbl->cj[found].excl &= ~(1U<<((inner_i<<na_cj_2log) + inner_e));
                    }
                }
            }
        }
    }
}

/* Add a new i-entry to the FEP list and copy the i-properties */
static gmx_inline void fep_list_new_nri_copy(t_nblist *nlist)
{
    /* Add a new i-entry */
    nlist->nri++;

    assert(nlist->nri < nlist->maxnri);

    /* Duplicate the last i-entry, except for jindex, which continues */
    nlist->iinr[nlist->nri]   = nlist->iinr[nlist->nri-1];
    nlist->shift[nlist->nri]  = nlist->shift[nlist->nri-1];
    nlist->gid[nlist->nri]    = nlist->gid[nlist->nri-1];
    nlist->jindex[nlist->nri] = nlist->nrj;
}

/* For load balancing of the free-energy lists over threads, we set
 * the maximum nrj size of an i-entry to 40. This leads to good
 * load balancing in the worst case scenario of a single perturbed
 * particle on 16 threads, while not introducing significant overhead.
 * Note that half of the perturbed pairs will anyhow end up in very small lists,
 * since non perturbed i-particles will see few perturbed j-particles).
 */
const int max_nrj_fep = 40;

/* Exclude the perturbed pairs from the Verlet list. This is only done to avoid
 * singularities for overlapping particles (0/0), since the charges and
 * LJ parameters have been zeroed in the nbnxn data structure.
 * Simultaneously make a group pair list for the perturbed pairs.
 */
static void make_fep_list(const nbnxn_search_t    nbs,
                          const nbnxn_atomdata_t *nbat,
                          nbnxn_pairlist_t       *nbl,
                          gmx_bool                bDiagRemoved,
                          nbnxn_ci_t             *nbl_ci,
                          const nbnxn_grid_t     *gridi,
                          const nbnxn_grid_t     *gridj,
                          t_nblist               *nlist)
{
    int      ci, cj_ind_start, cj_ind_end, cj_ind, cja, cjr;
    int      nri_max;
    int      ngid, gid_i = 0, gid_j, gid;
    int      egp_shift, egp_mask;
    int      gid_cj = 0;
    int      i, j, ind_i, ind_j, ai, aj;
    int      nri;
    gmx_bool bFEP_i, bFEP_i_all;

    if (nbl_ci->cj_ind_end == nbl_ci->cj_ind_start)
    {
        /* Empty list */
        return;
    }

    ci = nbl_ci->ci;

    cj_ind_start = nbl_ci->cj_ind_start;
    cj_ind_end   = nbl_ci->cj_ind_end;

    /* In worst case we have alternating energy groups
     * and create #atom-pair lists, which means we need the size
     * of a cluster pair (na_ci*na_cj) times the number of cj's.
     */
    nri_max = nbl->na_ci*nbl->na_cj*(cj_ind_end - cj_ind_start);
    if (nlist->nri + nri_max > nlist->maxnri)
    {
        nlist->maxnri = over_alloc_large(nlist->nri + nri_max);
        reallocate_nblist(nlist);
    }

    ngid = nbat->nenergrp;

    if (ngid*gridj->na_cj > sizeof(gid_cj)*8)
    {
        gmx_fatal(FARGS, "The Verlet scheme with %dx%d kernels and free-energy only supports up to %d energy groups",
                  gridi->na_c, gridj->na_cj, (sizeof(gid_cj)*8)/gridj->na_cj);
    }

    egp_shift = nbat->neg_2log;
    egp_mask  = (1<<nbat->neg_2log) - 1;

    /* Loop over the atoms in the i sub-cell */
    bFEP_i_all = TRUE;
    for (i = 0; i < nbl->na_ci; i++)
    {
        ind_i = ci*nbl->na_ci + i;
        ai    = nbs->a[ind_i];
        if (ai >= 0)
        {
            nri                  = nlist->nri;
            nlist->jindex[nri+1] = nlist->jindex[nri];
            nlist->iinr[nri]     = ai;
            /* The actual energy group pair index is set later */
            nlist->gid[nri]      = 0;
            nlist->shift[nri]    = nbl_ci->shift & NBNXN_CI_SHIFT;

            bFEP_i = gridi->fep[ci - gridi->cell0] & (1 << i);

            bFEP_i_all = bFEP_i_all && bFEP_i;

            if ((nlist->nrj + cj_ind_end - cj_ind_start)*nbl->na_cj > nlist->maxnrj)
            {
                nlist->maxnrj = over_alloc_small((nlist->nrj + cj_ind_end - cj_ind_start)*nbl->na_cj);
                srenew(nlist->jjnr,     nlist->maxnrj);
                srenew(nlist->excl_fep, nlist->maxnrj);
            }

            if (ngid > 1)
            {
                gid_i = (nbat->energrp[ci] >> (egp_shift*i)) & egp_mask;
            }

            for (cj_ind = cj_ind_start; cj_ind < cj_ind_end; cj_ind++)
            {
                unsigned int fep_cj;

                cja = nbl->cj[cj_ind].cj;

                if (gridj->na_cj == gridj->na_c)
                {
                    cjr    = cja - gridj->cell0;
                    fep_cj = gridj->fep[cjr];
                    if (ngid > 1)
                    {
                        gid_cj = nbat->energrp[cja];
                    }
                }
                else if (2*gridj->na_cj == gridj->na_c)
                {
                    cjr    = cja - gridj->cell0*2;
                    /* Extract half of the ci fep/energrp mask */
                    fep_cj = (gridj->fep[cjr>>1] >> ((cjr&1)*gridj->na_cj)) & ((1<<gridj->na_cj) - 1);
                    if (ngid > 1)
                    {
                        gid_cj = nbat->energrp[cja>>1] >> ((cja&1)*gridj->na_cj*egp_shift) & ((1<<(gridj->na_cj*egp_shift)) - 1);
                    }
                }
                else
                {
                    cjr    = cja - (gridj->cell0>>1);
                    /* Combine two ci fep masks/energrp */
                    fep_cj = gridj->fep[cjr*2] + (gridj->fep[cjr*2+1] << gridj->na_c);
                    if (ngid > 1)
                    {
                        gid_cj = nbat->energrp[cja*2] + (nbat->energrp[cja*2+1] << (gridj->na_c*egp_shift));
                    }
                }

                if (bFEP_i || fep_cj != 0)
                {
                    for (j = 0; j < nbl->na_cj; j++)
                    {
                        /* Is this interaction perturbed and not excluded? */
                        ind_j = cja*nbl->na_cj + j;
                        aj    = nbs->a[ind_j];
                        if (aj >= 0 &&
                            (bFEP_i || (fep_cj & (1 << j))) &&
                            (!bDiagRemoved || ind_j >= ind_i))
                        {
                            if (ngid > 1)
                            {
                                gid_j = (gid_cj >> (j*egp_shift)) & egp_mask;
                                gid   = GID(gid_i, gid_j, ngid);

                                if (nlist->nrj > nlist->jindex[nri] &&
                                    nlist->gid[nri] != gid)
                                {
                                    /* Energy group pair changed: new list */
                                    fep_list_new_nri_copy(nlist);
                                    nri = nlist->nri;
                                }
                                nlist->gid[nri] = gid;
                            }

                            if (nlist->nrj - nlist->jindex[nri] >= max_nrj_fep)
                            {
                                fep_list_new_nri_copy(nlist);
                                nri = nlist->nri;
                            }

                            /* Add it to the FEP list */
                            nlist->jjnr[nlist->nrj]     = aj;
                            nlist->excl_fep[nlist->nrj] = (nbl->cj[cj_ind].excl >> (i*nbl->na_cj + j)) & 1;
                            nlist->nrj++;

                            /* Exclude it from the normal list.
                             * Note that the charge has been set to zero,
                             * but we need to avoid 0/0, as perturbed atoms
                             * can be on top of each other.
                             */
                            nbl->cj[cj_ind].excl &= ~(1U << (i*nbl->na_cj + j));
                        }
                    }
                }
            }

            if (nlist->nrj > nlist->jindex[nri])
            {
                /* Actually add this new, non-empty, list */
                nlist->nri++;
                nlist->jindex[nlist->nri] = nlist->nrj;
            }
        }
    }

    if (bFEP_i_all)
    {
        /* All interactions are perturbed, we can skip this entry */
        nbl_ci->cj_ind_end = cj_ind_start;
    }
}

/* Return the index of atom a within a cluster */
static gmx_inline int cj_mod_cj4(int cj)
{
    return cj & (NBNXN_GPU_JGROUP_SIZE - 1);
}

/* Convert a j-cluster to a cj4 group */
static gmx_inline int cj_to_cj4(int cj)
{
    return cj >> NBNXN_GPU_JGROUP_SIZE_2LOG;
}

/* Return the index of an j-atom within a warp */
static gmx_inline int a_mod_wj(int a)
{
    return a & (NBNXN_GPU_CLUSTER_SIZE/2 - 1);
}

/* As make_fep_list above, but for super/sub lists. */
static void make_fep_list_supersub(const nbnxn_search_t    nbs,
                                   const nbnxn_atomdata_t *nbat,
                                   nbnxn_pairlist_t       *nbl,
                                   gmx_bool                bDiagRemoved,
                                   const nbnxn_sci_t      *nbl_sci,
                                   real                    shx,
                                   real                    shy,
                                   real                    shz,
                                   real                    rlist_fep2,
                                   const nbnxn_grid_t     *gridi,
                                   const nbnxn_grid_t     *gridj,
                                   t_nblist               *nlist)
{
    int                sci, cj4_ind_start, cj4_ind_end, cj4_ind, gcj, cjr;
    int                nri_max;
    int                c, c_abs;
    int                i, j, ind_i, ind_j, ai, aj;
    int                nri;
    gmx_bool           bFEP_i;
    real               xi, yi, zi;
    const nbnxn_cj4_t *cj4;

    if (nbl_sci->cj4_ind_end == nbl_sci->cj4_ind_start)
    {
        /* Empty list */
        return;
    }

    sci = nbl_sci->sci;

    cj4_ind_start = nbl_sci->cj4_ind_start;
    cj4_ind_end   = nbl_sci->cj4_ind_end;

    /* Here we process one super-cell, max #atoms na_sc, versus a list
     * cj4 entries, each with max NBNXN_GPU_JGROUP_SIZE cj's, each
     * of size na_cj atoms.
     * On the GPU we don't support energy groups (yet).
     * So for each of the na_sc i-atoms, we need max one FEP list
     * for each max_nrj_fep j-atoms.
     */
    nri_max = nbl->na_sc*nbl->na_cj*(1 + ((cj4_ind_end - cj4_ind_start)*NBNXN_GPU_JGROUP_SIZE)/max_nrj_fep);
    if (nlist->nri + nri_max > nlist->maxnri)
    {
        nlist->maxnri = over_alloc_large(nlist->nri + nri_max);
        reallocate_nblist(nlist);
    }

    /* Loop over the atoms in the i super-cluster */
    for (c = 0; c < GPU_NSUBCELL; c++)
    {
        c_abs = sci*GPU_NSUBCELL + c;

        for (i = 0; i < nbl->na_ci; i++)
        {
            ind_i = c_abs*nbl->na_ci + i;
            ai    = nbs->a[ind_i];
            if (ai >= 0)
            {
                nri                  = nlist->nri;
                nlist->jindex[nri+1] = nlist->jindex[nri];
                nlist->iinr[nri]     = ai;
                /* With GPUs, energy groups are not supported */
                nlist->gid[nri]      = 0;
                nlist->shift[nri]    = nbl_sci->shift & NBNXN_CI_SHIFT;

                bFEP_i = (gridi->fep[c_abs - gridi->cell0*GPU_NSUBCELL] & (1 << i));

                xi = nbat->x[ind_i*nbat->xstride+XX] + shx;
                yi = nbat->x[ind_i*nbat->xstride+YY] + shy;
                zi = nbat->x[ind_i*nbat->xstride+ZZ] + shz;

                if ((nlist->nrj + cj4_ind_end - cj4_ind_start)*NBNXN_GPU_JGROUP_SIZE*nbl->na_cj > nlist->maxnrj)
                {
                    nlist->maxnrj = over_alloc_small((nlist->nrj + cj4_ind_end - cj4_ind_start)*NBNXN_GPU_JGROUP_SIZE*nbl->na_cj);
                    srenew(nlist->jjnr,     nlist->maxnrj);
                    srenew(nlist->excl_fep, nlist->maxnrj);
                }

                for (cj4_ind = cj4_ind_start; cj4_ind < cj4_ind_end; cj4_ind++)
                {
                    cj4 = &nbl->cj4[cj4_ind];

                    for (gcj = 0; gcj < NBNXN_GPU_JGROUP_SIZE; gcj++)
                    {
                        unsigned int fep_cj;

                        if ((cj4->imei[0].imask & (1U << (gcj*GPU_NSUBCELL + c))) == 0)
                        {
                            /* Skip this ci for this cj */
                            continue;
                        }

                        cjr = cj4->cj[gcj] - gridj->cell0*GPU_NSUBCELL;

                        fep_cj = gridj->fep[cjr];

                        if (bFEP_i || fep_cj != 0)
                        {
                            for (j = 0; j < nbl->na_cj; j++)
                            {
                                /* Is this interaction perturbed and not excluded? */
                                ind_j = (gridj->cell0*GPU_NSUBCELL + cjr)*nbl->na_cj + j;
                                aj    = nbs->a[ind_j];
                                if (aj >= 0 &&
                                    (bFEP_i || (fep_cj & (1 << j))) &&
                                    (!bDiagRemoved || ind_j >= ind_i))
                                {
                                    nbnxn_excl_t *excl;
                                    int           excl_pair;
                                    unsigned int  excl_bit;
                                    real          dx, dy, dz;

                                    get_nbl_exclusions_1(nbl, cj4_ind, j>>2, &excl);

                                    excl_pair = a_mod_wj(j)*nbl->na_ci + i;
                                    excl_bit  = (1U << (gcj*GPU_NSUBCELL + c));

                                    dx = nbat->x[ind_j*nbat->xstride+XX] - xi;
                                    dy = nbat->x[ind_j*nbat->xstride+YY] - yi;
                                    dz = nbat->x[ind_j*nbat->xstride+ZZ] - zi;

                                    /* The unpruned GPU list has more than 2/3
                                     * of the atom pairs beyond rlist. Using
                                     * this list will cause a lot of overhead
                                     * in the CPU FEP kernels, especially
                                     * relative to the fast GPU kernels.
                                     * So we prune the FEP list here.
                                     */
                                    if (dx*dx + dy*dy + dz*dz < rlist_fep2)
                                    {
                                        if (nlist->nrj - nlist->jindex[nri] >= max_nrj_fep)
                                        {
                                            fep_list_new_nri_copy(nlist);
                                            nri = nlist->nri;
                                        }

                                        /* Add it to the FEP list */
                                        nlist->jjnr[nlist->nrj]     = aj;
                                        nlist->excl_fep[nlist->nrj] = (excl->pair[excl_pair] & excl_bit) ? 1 : 0;
                                        nlist->nrj++;
                                    }

                                    /* Exclude it from the normal list.
                                     * Note that the charge and LJ parameters have
                                     * been set to zero, but we need to avoid 0/0,
                                     * as perturbed atoms can be on top of each other.
                                     */
                                    excl->pair[excl_pair] &= ~excl_bit;
                                }
                            }

                            /* Note that we could mask out this pair in imask
                             * if all i- and/or all j-particles are perturbed.
                             * But since the perturbed pairs on the CPU will
                             * take an order of magnitude more time, the GPU
                             * will finish before the CPU and there is no gain.
                             */
                        }
                    }
                }

                if (nlist->nrj > nlist->jindex[nri])
                {
                    /* Actually add this new, non-empty, list */
                    nlist->nri++;
                    nlist->jindex[nlist->nri] = nlist->nrj;
                }
            }
        }
    }
}

/* Set all atom-pair exclusions from the topology stored in excl
 * as masks in the pair-list for i-super-cell entry nbl_sci
 */
static void set_sci_top_excls(const nbnxn_search_t nbs,
                              nbnxn_pairlist_t    *nbl,
                              gmx_bool             diagRemoved,
                              int                  na_c_2log,
                              const nbnxn_sci_t   *nbl_sci,
                              const t_blocka      *excl)
{
    const int    *cell;
    int           na_c;
    int           sci;
    int           cj_ind_first, cj_ind_last;
    int           cj_first, cj_last;
    int           ndirect;
    int           i, ai, aj, si, eind, ge, se;
    int           found, cj_ind_0, cj_ind_1, cj_ind_m;
    int           cj_m;
    gmx_bool      Found_si;
    int           si_ind;
    nbnxn_excl_t *nbl_excl;
    int           inner_i, inner_e, w;

    cell = nbs->cell;

    na_c = nbl->na_ci;

    if (nbl_sci->cj4_ind_end == nbl_sci->cj4_ind_start)
    {
        /* Empty list */
        return;
    }

    sci = nbl_sci->sci;

    cj_ind_first = nbl_sci->cj4_ind_start*NBNXN_GPU_JGROUP_SIZE;
    cj_ind_last  = nbl->work->cj_ind - 1;

    cj_first = nbl->cj4[nbl_sci->cj4_ind_start].cj[0];
    cj_last  = nbl_cj(nbl, cj_ind_last);

    /* Determine how many contiguous j-clusters we have starting
     * from the first i-cluster. This number can be used to directly
     * calculate j-cluster indices for excluded atoms.
     */
    ndirect = 0;
    while (cj_ind_first + ndirect <= cj_ind_last &&
           nbl_cj(nbl, cj_ind_first+ndirect) == sci*GPU_NSUBCELL + ndirect)
    {
        ndirect++;
    }

    /* Loop over the atoms in the i super-cell */
    for (i = 0; i < nbl->na_sc; i++)
    {
        ai = nbs->a[sci*nbl->na_sc+i];
        if (ai >= 0)
        {
            si  = (i>>na_c_2log);

            /* Loop over the topology-based exclusions for this i-atom */
            for (eind = excl->index[ai]; eind < excl->index[ai+1]; eind++)
            {
                aj = excl->a[eind];

                if (aj == ai)
                {
                    /* The self exclusion are already set, save some time */
                    continue;
                }

                ge = cell[aj];

                /* Without shifts we only calculate interactions j>i
                 * for one-way pair-lists.
                 */
                if (diagRemoved && ge <= sci*nbl->na_sc + i)
                {
                    continue;
                }

                se = ge>>na_c_2log;
                /* Could the cluster se be in our list? */
                if (se >= cj_first && se <= cj_last)
                {
                    if (se < cj_first + ndirect)
                    {
                        /* We can calculate cj_ind directly from se */
                        found = cj_ind_first + se - cj_first;
                    }
                    else
                    {
                        /* Search for se using bisection */
                        found    = -1;
                        cj_ind_0 = cj_ind_first + ndirect;
                        cj_ind_1 = cj_ind_last + 1;
                        while (found == -1 && cj_ind_0 < cj_ind_1)
                        {
                            cj_ind_m = (cj_ind_0 + cj_ind_1)>>1;

                            cj_m = nbl_cj(nbl, cj_ind_m);

                            if (se == cj_m)
                            {
                                found = cj_ind_m;
                            }
                            else if (se < cj_m)
                            {
                                cj_ind_1 = cj_ind_m;
                            }
                            else
                            {
                                cj_ind_0 = cj_ind_m + 1;
                            }
                        }
                    }

                    if (found >= 0)
                    {
                        inner_i = i  - si*na_c;
                        inner_e = ge - se*na_c;

                        if (nbl_imask0(nbl, found) & (1U << (cj_mod_cj4(found)*GPU_NSUBCELL + si)))
                        {
                            w       = (inner_e >> 2);

                            get_nbl_exclusions_1(nbl, cj_to_cj4(found), w, &nbl_excl);

                            nbl_excl->pair[a_mod_wj(inner_e)*nbl->na_ci+inner_i] &=
                                ~(1U << (cj_mod_cj4(found)*GPU_NSUBCELL + si));
                        }
                    }
                }
            }
        }
    }
}

/* Reallocate the simple ci list for at least n entries */
static void nb_realloc_ci(nbnxn_pairlist_t *nbl, int n)
{
    nbl->ci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->ci,
                       nbl->nci*sizeof(*nbl->ci),
                       nbl->ci_nalloc*sizeof(*nbl->ci),
                       nbl->alloc, nbl->free);
}

/* Reallocate the super-cell sci list for at least n entries */
static void nb_realloc_sci(nbnxn_pairlist_t *nbl, int n)
{
    nbl->sci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->sci,
                       nbl->nsci*sizeof(*nbl->sci),
                       nbl->sci_nalloc*sizeof(*nbl->sci),
                       nbl->alloc, nbl->free);
}

/* Make a new ci entry at index nbl->nci */
static void new_ci_entry(nbnxn_pairlist_t *nbl, int ci, int shift, int flags)
{
    if (nbl->nci + 1 > nbl->ci_nalloc)
    {
        nb_realloc_ci(nbl, nbl->nci+1);
    }
    nbl->ci[nbl->nci].ci            = ci;
    nbl->ci[nbl->nci].shift         = shift;
    /* Store the interaction flags along with the shift */
    nbl->ci[nbl->nci].shift        |= flags;
    nbl->ci[nbl->nci].cj_ind_start  = nbl->ncj;
    nbl->ci[nbl->nci].cj_ind_end    = nbl->ncj;
}

/* Make a new sci entry at index nbl->nsci */
static void new_sci_entry(nbnxn_pairlist_t *nbl, int sci, int shift)
{
    if (nbl->nsci + 1 > nbl->sci_nalloc)
    {
        nb_realloc_sci(nbl, nbl->nsci+1);
    }
    nbl->sci[nbl->nsci].sci           = sci;
    nbl->sci[nbl->nsci].shift         = shift;
    nbl->sci[nbl->nsci].cj4_ind_start = nbl->ncj4;
    nbl->sci[nbl->nsci].cj4_ind_end   = nbl->ncj4;
}

/* Sort the simple j-list cj on exclusions.
 * Entries with exclusions will all be sorted to the beginning of the list.
 */
static void sort_cj_excl(nbnxn_cj_t *cj, int ncj,
                         nbnxn_list_work_t *work)
{
    int jnew, j;

    if (ncj > work->cj_nalloc)
    {
        work->cj_nalloc = over_alloc_large(ncj);
        srenew(work->cj, work->cj_nalloc);
    }

    /* Make a list of the j-cells involving exclusions */
    jnew = 0;
    for (j = 0; j < ncj; j++)
    {
        if (cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            work->cj[jnew++] = cj[j];
        }
    }
    /* Check if there are exclusions at all or not just the first entry */
    if (!((jnew == 0) ||
          (jnew == 1 && cj[0].excl != NBNXN_INTERACTION_MASK_ALL)))
    {
        for (j = 0; j < ncj; j++)
        {
            if (cj[j].excl == NBNXN_INTERACTION_MASK_ALL)
            {
                work->cj[jnew++] = cj[j];
            }
        }
        for (j = 0; j < ncj; j++)
        {
            cj[j] = work->cj[j];
        }
    }
}

/* Close this simple list i entry */
static void close_ci_entry_simple(nbnxn_pairlist_t *nbl)
{
    int jlen;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    jlen = nbl->ci[nbl->nci].cj_ind_end - nbl->ci[nbl->nci].cj_ind_start;
    if (jlen > 0)
    {
        sort_cj_excl(nbl->cj+nbl->ci[nbl->nci].cj_ind_start, jlen, nbl->work);

        /* The counts below are used for non-bonded pair/flop counts
         * and should therefore match the available kernel setups.
         */
        if (!(nbl->ci[nbl->nci].shift & NBNXN_CI_DO_COUL(0)))
        {
            nbl->work->ncj_noq += jlen;
        }
        else if ((nbl->ci[nbl->nci].shift & NBNXN_CI_HALF_LJ(0)) ||
                 !(nbl->ci[nbl->nci].shift & NBNXN_CI_DO_LJ(0)))
        {
            nbl->work->ncj_hlj += jlen;
        }

        nbl->nci++;
    }
}

/* Split sci entry for load balancing on the GPU.
 * Splitting ensures we have enough lists to fully utilize the whole GPU.
 * With progBal we generate progressively smaller lists, which improves
 * load balancing. As we only know the current count on our own thread,
 * we will need to estimate the current total amount of i-entries.
 * As the lists get concatenated later, this estimate depends
 * both on nthread and our own thread index.
 */
static void split_sci_entry(nbnxn_pairlist_t *nbl,
                            int nsp_target_av,
                            gmx_bool progBal, int nsp_tot_est,
                            int thread, int nthread)
{
    int nsp_est;
    int nsp_max;
    int cj4_start, cj4_end, j4len, cj4;
    int sci;
    int nsp, nsp_sci, nsp_cj4, nsp_cj4_e, nsp_cj4_p;
    int p;

    if (progBal)
    {
        /* Estimate the total numbers of ci's of the nblist combined
         * over all threads using the target number of ci's.
         */
        nsp_est = (nsp_tot_est*thread)/nthread + nbl->nci_tot;

        /* The first ci blocks should be larger, to avoid overhead.
         * The last ci blocks should be smaller, to improve load balancing.
         * The factor 3/2 makes the first block 3/2 times the target average
         * and ensures that the total number of blocks end up equal to
         * that with of equally sized blocks of size nsp_target_av.
         */
        nsp_max = nsp_target_av*nsp_tot_est*3/(2*(nsp_est + nsp_tot_est));
    }
    else
    {
        nsp_max = nsp_target_av;
    }

    /* Since nsp_max is a maximum/cut-off (this avoids high outliers,
     * which lead to load imbalance), not an average, we add half the
     * number of pairs in a cj4 block to get the average about right.
     */
    nsp_max += GPU_NSUBCELL*NBNXN_GPU_JGROUP_SIZE/2;

    cj4_start = nbl->sci[nbl->nsci-1].cj4_ind_start;
    cj4_end   = nbl->sci[nbl->nsci-1].cj4_ind_end;
    j4len     = cj4_end - cj4_start;

    if (j4len > 1 && j4len*GPU_NSUBCELL*NBNXN_GPU_JGROUP_SIZE > nsp_max)
    {
        /* Remove the last ci entry and process the cj4's again */
        nbl->nsci -= 1;

        sci        = nbl->nsci;
        nsp        = 0;
        nsp_sci    = 0;
        nsp_cj4_e  = 0;
        nsp_cj4    = 0;
        for (cj4 = cj4_start; cj4 < cj4_end; cj4++)
        {
            nsp_cj4_p = nsp_cj4;
            /* Count the number of cluster pairs in this cj4 group */
            nsp_cj4   = 0;
            for (p = 0; p < GPU_NSUBCELL*NBNXN_GPU_JGROUP_SIZE; p++)
            {
                nsp_cj4 += (nbl->cj4[cj4].imei[0].imask >> p) & 1;
            }

            /* Check if we should split at this cj4 to get a list of size nsp */
            if (nsp > 0 && nsp + nsp_cj4 > nsp_max)
            {
                /* Split the list at cj4 */
                nbl->sci[sci].cj4_ind_end = cj4;
                /* Create a new sci entry */
                sci++;
                nbl->nsci++;
                if (nbl->nsci+1 > nbl->sci_nalloc)
                {
                    nb_realloc_sci(nbl, nbl->nsci+1);
                }
                nbl->sci[sci].sci           = nbl->sci[nbl->nsci-1].sci;
                nbl->sci[sci].shift         = nbl->sci[nbl->nsci-1].shift;
                nbl->sci[sci].cj4_ind_start = cj4;
                nsp_sci                     = nsp;
                nsp_cj4_e                   = nsp_cj4_p;
                nsp                         = 0;
            }
            nsp += nsp_cj4;
        }

        /* Put the remaining cj4's in the last sci entry */
        nbl->sci[sci].cj4_ind_end = cj4_end;

        /* Possibly balance out the last two sci's
         * by moving the last cj4 of the second last sci.
         */
        if (nsp_sci - nsp_cj4_e >= nsp + nsp_cj4_e)
        {
            nbl->sci[sci-1].cj4_ind_end--;
            nbl->sci[sci].cj4_ind_start--;
        }

        nbl->nsci++;
    }
}

/* Clost this super/sub list i entry */
static void close_ci_entry_supersub(nbnxn_pairlist_t *nbl,
                                    int nsp_max_av,
                                    gmx_bool progBal, int nsp_tot_est,
                                    int thread, int nthread)
{
    int j4len, tlen;
    int nb, b;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    j4len = nbl->sci[nbl->nsci].cj4_ind_end - nbl->sci[nbl->nsci].cj4_ind_start;
    if (j4len > 0)
    {
        /* We can only have complete blocks of 4 j-entries in a list,
         * so round the count up before closing.
         */
        nbl->ncj4         = ((nbl->work->cj_ind + NBNXN_GPU_JGROUP_SIZE - 1) >> NBNXN_GPU_JGROUP_SIZE_2LOG);
        nbl->work->cj_ind = nbl->ncj4*NBNXN_GPU_JGROUP_SIZE;

        nbl->nsci++;

        if (nsp_max_av > 0)
        {
            /* Measure the size of the new entry and potentially split it */
            split_sci_entry(nbl, nsp_max_av, progBal, nsp_tot_est,
                            thread, nthread);
        }
    }
}

/* Syncs the working array before adding another grid pair to the list */
static void sync_work(nbnxn_pairlist_t *nbl)
{
    if (!nbl->bSimple)
    {
        nbl->work->cj_ind   = nbl->ncj4*NBNXN_GPU_JGROUP_SIZE;
        nbl->work->cj4_init = nbl->ncj4;
    }
}

/* Clears an nbnxn_pairlist_t data structure */
static void clear_pairlist(nbnxn_pairlist_t *nbl)
{
    nbl->nci           = 0;
    nbl->nsci          = 0;
    nbl->ncj           = 0;
    nbl->ncj4          = 0;
    nbl->nci_tot       = 0;
    nbl->nexcl         = 1;

    nbl->work->ncj_noq = 0;
    nbl->work->ncj_hlj = 0;
}

/* Clears a group scheme pair list */
static void clear_pairlist_fep(t_nblist *nl)
{
    nl->nri = 0;
    nl->nrj = 0;
    if (nl->jindex == NULL)
    {
        snew(nl->jindex, 1);
    }
    nl->jindex[0] = 0;
}

/* Sets a simple list i-cell bounding box, including PBC shift */
static gmx_inline void set_icell_bb_simple(const nbnxn_bb_t *bb, int ci,
                                           real shx, real shy, real shz,
                                           nbnxn_bb_t *bb_ci)
{
    bb_ci->lower[BB_X] = bb[ci].lower[BB_X] + shx;
    bb_ci->lower[BB_Y] = bb[ci].lower[BB_Y] + shy;
    bb_ci->lower[BB_Z] = bb[ci].lower[BB_Z] + shz;
    bb_ci->upper[BB_X] = bb[ci].upper[BB_X] + shx;
    bb_ci->upper[BB_Y] = bb[ci].upper[BB_Y] + shy;
    bb_ci->upper[BB_Z] = bb[ci].upper[BB_Z] + shz;
}

#ifdef NBNXN_BBXXXX
/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
static void set_icell_bbxxxx_supersub(const float *bb, int ci,
                                      real shx, real shy, real shz,
                                      float *bb_ci)
{
    int ia, m, i;

    ia = ci*(GPU_NSUBCELL>>STRIDE_PBB_2LOG)*NNBSBB_XXXX;
    for (m = 0; m < (GPU_NSUBCELL>>STRIDE_PBB_2LOG)*NNBSBB_XXXX; m += NNBSBB_XXXX)
    {
        for (i = 0; i < STRIDE_PBB; i++)
        {
            bb_ci[m+0*STRIDE_PBB+i] = bb[ia+m+0*STRIDE_PBB+i] + shx;
            bb_ci[m+1*STRIDE_PBB+i] = bb[ia+m+1*STRIDE_PBB+i] + shy;
            bb_ci[m+2*STRIDE_PBB+i] = bb[ia+m+2*STRIDE_PBB+i] + shz;
            bb_ci[m+3*STRIDE_PBB+i] = bb[ia+m+3*STRIDE_PBB+i] + shx;
            bb_ci[m+4*STRIDE_PBB+i] = bb[ia+m+4*STRIDE_PBB+i] + shy;
            bb_ci[m+5*STRIDE_PBB+i] = bb[ia+m+5*STRIDE_PBB+i] + shz;
        }
    }
}
#endif

/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
static void set_icell_bb_supersub(const nbnxn_bb_t *bb, int ci,
                                  real shx, real shy, real shz,
                                  nbnxn_bb_t *bb_ci)
{
    int i;

    for (i = 0; i < GPU_NSUBCELL; i++)
    {
        set_icell_bb_simple(bb, ci*GPU_NSUBCELL+i,
                            shx, shy, shz,
                            &bb_ci[i]);
    }
}

/* Copies PBC shifted i-cell atom coordinates x,y,z to working array */
static void icell_set_x_simple(int ci,
                               real shx, real shy, real shz,
                               int gmx_unused na_c,
                               int stride, const real *x,
                               nbnxn_list_work_t *work)
{
    int  ia, i;

    ia = ci*NBNXN_CPU_CLUSTER_I_SIZE;

    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        work->x_ci[i*STRIDE_XYZ+XX] = x[(ia+i)*stride+XX] + shx;
        work->x_ci[i*STRIDE_XYZ+YY] = x[(ia+i)*stride+YY] + shy;
        work->x_ci[i*STRIDE_XYZ+ZZ] = x[(ia+i)*stride+ZZ] + shz;
    }
}

/* Copies PBC shifted super-cell atom coordinates x,y,z to working array */
static void icell_set_x_supersub(int ci,
                                 real shx, real shy, real shz,
                                 int na_c,
                                 int stride, const real *x,
                                 nbnxn_list_work_t *work)
{
    int   ia, i;
    real *x_ci;

    x_ci = work->x_ci;

    ia = ci*GPU_NSUBCELL*na_c;
    for (i = 0; i < GPU_NSUBCELL*na_c; i++)
    {
        x_ci[i*DIM + XX] = x[(ia+i)*stride + XX] + shx;
        x_ci[i*DIM + YY] = x[(ia+i)*stride + YY] + shy;
        x_ci[i*DIM + ZZ] = x[(ia+i)*stride + ZZ] + shz;
    }
}

#ifdef NBNXN_SEARCH_BB_SIMD4
/* Copies PBC shifted super-cell packed atom coordinates to working array */
static void icell_set_x_supersub_simd4(int ci,
                                       real shx, real shy, real shz,
                                       int na_c,
                                       int stride, const real *x,
                                       nbnxn_list_work_t *work)
{
    int   si, io, ia, i, j;
    real *x_ci;

    x_ci = work->x_ci;

    for (si = 0; si < GPU_NSUBCELL; si++)
    {
        for (i = 0; i < na_c; i += STRIDE_PBB)
        {
            io = si*na_c + i;
            ia = ci*GPU_NSUBCELL*na_c + io;
            for (j = 0; j < STRIDE_PBB; j++)
            {
                x_ci[io*DIM + j + XX*STRIDE_PBB] = x[(ia+j)*stride+XX] + shx;
                x_ci[io*DIM + j + YY*STRIDE_PBB] = x[(ia+j)*stride+YY] + shy;
                x_ci[io*DIM + j + ZZ*STRIDE_PBB] = x[(ia+j)*stride+ZZ] + shz;
            }
        }
    }
}
#endif

static real minimum_subgrid_size_xy(const nbnxn_grid_t *grid)
{
    if (grid->bSimple)
    {
        return min(grid->sx, grid->sy);
    }
    else
    {
        return min(grid->sx/GPU_NSUBCELL_X, grid->sy/GPU_NSUBCELL_Y);
    }
}

static real effective_buffer_1x1_vs_MxN(const nbnxn_grid_t *gridi,
                                        const nbnxn_grid_t *gridj)
{
    const real eff_1x1_buffer_fac_overest = 0.1;

    /* Determine an atom-pair list cut-off buffer size for atom pairs,
     * to be added to rlist (including buffer) used for MxN.
     * This is for converting an MxN list to a 1x1 list. This means we can't
     * use the normal buffer estimate, as we have an MxN list in which
     * some atom pairs beyond rlist are missing. We want to capture
     * the beneficial effect of buffering by extra pairs just outside rlist,
     * while removing the useless pairs that are further away from rlist.
     * (Also the buffer could have been set manually not using the estimate.)
     * This buffer size is an overestimate.
     * We add 10% of the smallest grid sub-cell dimensions.
     * Note that the z-size differs per cell and we don't use this,
     * so we overestimate.
     * With PME, the 10% value gives a buffer that is somewhat larger
     * than the effective buffer with a tolerance of 0.005 kJ/mol/ps.
     * Smaller tolerances or using RF lead to a smaller effective buffer,
     * so 10% gives a safe overestimate.
     */
    return eff_1x1_buffer_fac_overest*(minimum_subgrid_size_xy(gridi) +
                                       minimum_subgrid_size_xy(gridj));
}

/* Clusters at the cut-off only increase rlist by 60% of their size */
static real nbnxn_rlist_inc_outside_fac = 0.6;

/* Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 */
real nbnxn_get_rlist_effective_inc(int cluster_size_j, real atom_density)
{
    int  cluster_size_i;
    real vol_inc_i, vol_inc_j;

    /* We should get this from the setup, but currently it's the same for
     * all setups, including GPUs.
     */
    cluster_size_i = NBNXN_CPU_CLUSTER_I_SIZE;

    vol_inc_i = (cluster_size_i - 1)/atom_density;
    vol_inc_j = (cluster_size_j - 1)/atom_density;

    return nbnxn_rlist_inc_outside_fac*pow(vol_inc_i + vol_inc_j, 1.0/3.0);
}

/* Estimates the interaction volume^2 for non-local interactions */
static real nonlocal_vol2(const gmx_domdec_zones_t *zones, rvec ls, real r)
{
    int  z, d;
    real cl, ca, za;
    real vold_est;
    real vol2_est_tot;

    vol2_est_tot = 0;

    /* Here we simply add up the volumes of 1, 2 or 3 1D decomposition
     * not home interaction volume^2. As these volumes are not additive,
     * this is an overestimate, but it would only be significant in the limit
     * of small cells, where we anyhow need to split the lists into
     * as small parts as possible.
     */

    for (z = 0; z < zones->n; z++)
    {
        if (zones->shift[z][XX] + zones->shift[z][YY] + zones->shift[z][ZZ] == 1)
        {
            cl = 0;
            ca = 1;
            za = 1;
            for (d = 0; d < DIM; d++)
            {
                if (zones->shift[z][d] == 0)
                {
                    cl += 0.5*ls[d];
                    ca *= ls[d];
                    za *= zones->size[z].x1[d] - zones->size[z].x0[d];
                }
            }

            /* 4 octants of a sphere */
            vold_est  = 0.25*M_PI*r*r*r*r;
            /* 4 quarter pie slices on the edges */
            vold_est += 4*cl*M_PI/6.0*r*r*r;
            /* One rectangular volume on a face */
            vold_est += ca*0.5*r*r;

            vol2_est_tot += vold_est*za;
        }
    }

    return vol2_est_tot;
}

/* Estimates the average size of a full j-list for super/sub setup */
static void get_nsubpair_target(const nbnxn_search_t  nbs,
                                int                   iloc,
                                real                  rlist,
                                int                   min_ci_balanced,
                                int                  *nsubpair_target,
                                int                  *nsubpair_tot_est)
{
    /* The target value of 36 seems to be the optimum for Kepler.
     * Maxwell is less sensitive to the exact value.
     */
    const int           nsubpair_target_min = 36;
    const nbnxn_grid_t *grid;
    rvec                ls;
    real                xy_diag2, r_eff_sup, vol_est, nsp_est, nsp_est_nl;
    int                 nsubpair_max;

    grid = &nbs->grid[0];

    if (min_ci_balanced <= 0 || grid->nc >= min_ci_balanced || grid->nc == 0)
    {
        /* We don't need to balance the list sizes */
        *nsubpair_target  = 0;
        *nsubpair_tot_est = 0;

        return;
    }

    ls[XX] = (grid->c1[XX] - grid->c0[XX])/(grid->ncx*GPU_NSUBCELL_X);
    ls[YY] = (grid->c1[YY] - grid->c0[YY])/(grid->ncy*GPU_NSUBCELL_Y);
    ls[ZZ] = grid->na_c/(grid->atom_density*ls[XX]*ls[YY]);

    /* The average squared length of the diagonal of a sub cell */
    xy_diag2 = ls[XX]*ls[XX] + ls[YY]*ls[YY] + ls[ZZ]*ls[ZZ];

    /* The formulas below are a heuristic estimate of the average nsj per si*/
    r_eff_sup = rlist + nbnxn_rlist_inc_outside_fac*sqr((grid->na_c - 1.0)/grid->na_c)*sqrt(xy_diag2/3);

    if (!nbs->DomDec || nbs->zones->n == 1)
    {
        nsp_est_nl = 0;
    }
    else
    {
        nsp_est_nl =
            sqr(grid->atom_density/grid->na_c)*
            nonlocal_vol2(nbs->zones, ls, r_eff_sup);
    }

    if (LOCAL_I(iloc))
    {
        /* Sub-cell interacts with itself */
        vol_est  = ls[XX]*ls[YY]*ls[ZZ];
        /* 6/2 rectangular volume on the faces */
        vol_est += (ls[XX]*ls[YY] + ls[XX]*ls[ZZ] + ls[YY]*ls[ZZ])*r_eff_sup;
        /* 12/2 quarter pie slices on the edges */
        vol_est += 2*(ls[XX] + ls[YY] + ls[ZZ])*0.25*M_PI*sqr(r_eff_sup);
        /* 4 octants of a sphere */
        vol_est += 0.5*4.0/3.0*M_PI*pow(r_eff_sup, 3);

        /* Estimate the number of cluster pairs as the local number of
         * clusters times the volume they interact with times the density.
         */
        nsp_est = grid->nsubc_tot*vol_est*grid->atom_density/grid->na_c;

        /* Subtract the non-local pair count */
        nsp_est -= nsp_est_nl;

        /* For small cut-offs nsp_est will be an underesimate.
         * With DD nsp_est_nl is an overestimate so nsp_est can get negative.
         * So to avoid too small or negative nsp_est we set a minimum of
         * all cells interacting with all 3^3 direct neighbors (3^3-1)/2+1=14.
         * This might be a slight overestimate for small non-periodic groups of
         * atoms as will occur for a local domain with DD, but for small
         * groups of atoms we'll anyhow be limited by nsubpair_target_min,
         * so this overestimation will not matter.
         */
        nsp_est = max(nsp_est, grid->nsubc_tot*14.0);

        if (debug)
        {
            fprintf(debug, "nsp_est local %5.1f non-local %5.1f\n",
                    nsp_est, nsp_est_nl);
        }
    }
    else
    {
        nsp_est = nsp_est_nl;
    }

    /* Thus the (average) maximum j-list size should be as follows.
     * Since there is overhead, we shouldn't make the lists too small
     * (and we can't chop up j-groups) so we use a minimum target size of 36.
     */
    *nsubpair_target  = max(nsubpair_target_min,
                            (int)(nsp_est/min_ci_balanced + 0.5));
    *nsubpair_tot_est = (int)nsp_est;

    if (debug)
    {
        fprintf(debug, "nbl nsp estimate %.1f, nsubpair_target %d\n",
                nsp_est, *nsubpair_target);
    }
}

/* Debug list print function */
static void print_nblist_ci_cj(FILE *fp, const nbnxn_pairlist_t *nbl)
{
    int i, j;

    for (i = 0; i < nbl->nci; i++)
    {
        fprintf(fp, "ci %4d  shift %2d  ncj %3d\n",
                nbl->ci[i].ci, nbl->ci[i].shift,
                nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start);

        for (j = nbl->ci[i].cj_ind_start; j < nbl->ci[i].cj_ind_end; j++)
        {
            fprintf(fp, "  cj %5d  imask %x\n",
                    nbl->cj[j].cj,
                    nbl->cj[j].excl);
        }
    }
}

/* Debug list print function */
static void print_nblist_sci_cj(FILE *fp, const nbnxn_pairlist_t *nbl)
{
    int i, j4, j, ncp, si;

    for (i = 0; i < nbl->nsci; i++)
    {
        fprintf(fp, "ci %4d  shift %2d  ncj4 %2d\n",
                nbl->sci[i].sci, nbl->sci[i].shift,
                nbl->sci[i].cj4_ind_end - nbl->sci[i].cj4_ind_start);

        ncp = 0;
        for (j4 = nbl->sci[i].cj4_ind_start; j4 < nbl->sci[i].cj4_ind_end; j4++)
        {
            for (j = 0; j < NBNXN_GPU_JGROUP_SIZE; j++)
            {
                fprintf(fp, "  sj %5d  imask %x\n",
                        nbl->cj4[j4].cj[j],
                        nbl->cj4[j4].imei[0].imask);
                for (si = 0; si < GPU_NSUBCELL; si++)
                {
                    if (nbl->cj4[j4].imei[0].imask & (1U << (j*GPU_NSUBCELL + si)))
                    {
                        ncp++;
                    }
                }
            }
        }
        fprintf(fp, "ci %4d  shift %2d  ncj4 %2d ncp %3d\n",
                nbl->sci[i].sci, nbl->sci[i].shift,
                nbl->sci[i].cj4_ind_end - nbl->sci[i].cj4_ind_start,
                ncp);
    }
}

/* Combine pair lists *nbl generated on multiple threads nblc */
static void combine_nblists(int nnbl, nbnxn_pairlist_t **nbl,
                            nbnxn_pairlist_t *nblc)
{
    int nsci, ncj4, nexcl;
    int n, i;
    int nthreads gmx_unused;

    if (nblc->bSimple)
    {
        gmx_incons("combine_nblists does not support simple lists");
    }

    nsci  = nblc->nsci;
    ncj4  = nblc->ncj4;
    nexcl = nblc->nexcl;
    for (i = 0; i < nnbl; i++)
    {
        nsci  += nbl[i]->nsci;
        ncj4  += nbl[i]->ncj4;
        nexcl += nbl[i]->nexcl;
    }

    if (nsci > nblc->sci_nalloc)
    {
        nb_realloc_sci(nblc, nsci);
    }
    if (ncj4 > nblc->cj4_nalloc)
    {
        nblc->cj4_nalloc = over_alloc_small(ncj4);
        nbnxn_realloc_void((void **)&nblc->cj4,
                           nblc->ncj4*sizeof(*nblc->cj4),
                           nblc->cj4_nalloc*sizeof(*nblc->cj4),
                           nblc->alloc, nblc->free);
    }
    if (nexcl > nblc->excl_nalloc)
    {
        nblc->excl_nalloc = over_alloc_small(nexcl);
        nbnxn_realloc_void((void **)&nblc->excl,
                           nblc->nexcl*sizeof(*nblc->excl),
                           nblc->excl_nalloc*sizeof(*nblc->excl),
                           nblc->alloc, nblc->free);
    }

    /* Each thread should copy its own data to the combined arrays,
     * as otherwise data will go back and forth between different caches.
     */
    nthreads = gmx_omp_nthreads_get(emntPairsearch);
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (n = 0; n < nnbl; n++)
    {
        int                     sci_offset;
        int                     cj4_offset;
        int                     ci_offset;
        int                     excl_offset;
        int                     i, j4;
        const nbnxn_pairlist_t *nbli;

        /* Determine the offset in the combined data for our thread */
        sci_offset  = nblc->nsci;
        cj4_offset  = nblc->ncj4;
        ci_offset   = nblc->nci_tot;
        excl_offset = nblc->nexcl;

        for (i = 0; i < n; i++)
        {
            sci_offset  += nbl[i]->nsci;
            cj4_offset  += nbl[i]->ncj4;
            ci_offset   += nbl[i]->nci_tot;
            excl_offset += nbl[i]->nexcl;
        }

        nbli = nbl[n];

        for (i = 0; i < nbli->nsci; i++)
        {
            nblc->sci[sci_offset+i]                = nbli->sci[i];
            nblc->sci[sci_offset+i].cj4_ind_start += cj4_offset;
            nblc->sci[sci_offset+i].cj4_ind_end   += cj4_offset;
        }

        for (j4 = 0; j4 < nbli->ncj4; j4++)
        {
            nblc->cj4[cj4_offset+j4]                   = nbli->cj4[j4];
            nblc->cj4[cj4_offset+j4].imei[0].excl_ind += excl_offset;
            nblc->cj4[cj4_offset+j4].imei[1].excl_ind += excl_offset;
        }

        for (j4 = 0; j4 < nbli->nexcl; j4++)
        {
            nblc->excl[excl_offset+j4] = nbli->excl[j4];
        }
    }

    for (n = 0; n < nnbl; n++)
    {
        nblc->nsci    += nbl[n]->nsci;
        nblc->ncj4    += nbl[n]->ncj4;
        nblc->nci_tot += nbl[n]->nci_tot;
        nblc->nexcl   += nbl[n]->nexcl;
    }
}

static void balance_fep_lists(const nbnxn_search_t  nbs,
                              nbnxn_pairlist_set_t *nbl_lists)
{
    int       nnbl, th;
    int       nri_tot, nrj_tot, nrj_target;
    int       th_dest;
    t_nblist *nbld;

    nnbl = nbl_lists->nnbl;

    if (nnbl == 1)
    {
        /* Nothing to balance */
        return;
    }

    /* Count the total i-lists and pairs */
    nri_tot = 0;
    nrj_tot = 0;
    for (th = 0; th < nnbl; th++)
    {
        nri_tot += nbl_lists->nbl_fep[th]->nri;
        nrj_tot += nbl_lists->nbl_fep[th]->nrj;
    }

    nrj_target = (nrj_tot + nnbl - 1)/nnbl;

    assert(gmx_omp_nthreads_get(emntNonbonded) == nnbl);

#pragma omp parallel for schedule(static) num_threads(nnbl)
    for (th = 0; th < nnbl; th++)
    {
        t_nblist *nbl;

        nbl = nbs->work[th].nbl_fep;

        /* Note that here we allocate for the total size, instead of
         * a per-thread esimate (which is hard to obtain).
         */
        if (nri_tot > nbl->maxnri)
        {
            nbl->maxnri = over_alloc_large(nri_tot);
            reallocate_nblist(nbl);
        }
        if (nri_tot > nbl->maxnri || nrj_tot > nbl->maxnrj)
        {
            nbl->maxnrj = over_alloc_small(nrj_tot);
            srenew(nbl->jjnr, nbl->maxnrj);
            srenew(nbl->excl_fep, nbl->maxnrj);
        }

        clear_pairlist_fep(nbl);
    }

    /* Loop over the source lists and assign and copy i-entries */
    th_dest = 0;
    nbld    = nbs->work[th_dest].nbl_fep;
    for (th = 0; th < nnbl; th++)
    {
        t_nblist *nbls;
        int       i, j;

        nbls = nbl_lists->nbl_fep[th];

        for (i = 0; i < nbls->nri; i++)
        {
            int nrj;

            /* The number of pairs in this i-entry */
            nrj = nbls->jindex[i+1] - nbls->jindex[i];

            /* Decide if list th_dest is too large and we should procede
             * to the next destination list.
             */
            if (th_dest+1 < nnbl && nbld->nrj > 0 &&
                nbld->nrj + nrj - nrj_target > nrj_target - nbld->nrj)
            {
                th_dest++;
                nbld = nbs->work[th_dest].nbl_fep;
            }

            nbld->iinr[nbld->nri]  = nbls->iinr[i];
            nbld->gid[nbld->nri]   = nbls->gid[i];
            nbld->shift[nbld->nri] = nbls->shift[i];

            for (j = nbls->jindex[i]; j < nbls->jindex[i+1]; j++)
            {
                nbld->jjnr[nbld->nrj]     = nbls->jjnr[j];
                nbld->excl_fep[nbld->nrj] = nbls->excl_fep[j];
                nbld->nrj++;
            }
            nbld->nri++;
            nbld->jindex[nbld->nri] = nbld->nrj;
        }
    }

    /* Swap the list pointers */
    for (th = 0; th < nnbl; th++)
    {
        t_nblist *nbl_tmp;

        nbl_tmp                = nbl_lists->nbl_fep[th];
        nbl_lists->nbl_fep[th] = nbs->work[th].nbl_fep;
        nbs->work[th].nbl_fep  = nbl_tmp;

        if (debug)
        {
            fprintf(debug, "nbl_fep[%d] nri %4d nrj %4d\n",
                    th,
                    nbl_lists->nbl_fep[th]->nri,
                    nbl_lists->nbl_fep[th]->nrj);
        }
    }
}

/* Returns the next ci to be processes by our thread */
static gmx_bool next_ci(const nbnxn_grid_t *grid,
                        int conv,
                        int nth, int ci_block,
                        int *ci_x, int *ci_y,
                        int *ci_b, int *ci)
{
    (*ci_b)++;
    (*ci)++;

    if (*ci_b == ci_block)
    {
        /* Jump to the next block assigned to this task */
        *ci   += (nth - 1)*ci_block;
        *ci_b  = 0;
    }

    if (*ci >= grid->nc*conv)
    {
        return FALSE;
    }

    while (*ci >= grid->cxy_ind[*ci_x*grid->ncy + *ci_y + 1]*conv)
    {
        *ci_y += 1;
        if (*ci_y == grid->ncy)
        {
            *ci_x += 1;
            *ci_y  = 0;
        }
    }

    return TRUE;
}

/* Returns the distance^2 for which we put cell pairs in the list
 * without checking atom pair distances. This is usually < rlist^2.
 */
static float boundingbox_only_distance2(const nbnxn_grid_t *gridi,
                                        const nbnxn_grid_t *gridj,
                                        real                rlist,
                                        gmx_bool            simple)
{
    /* If the distance between two sub-cell bounding boxes is less
     * than this distance, do not check the distance between
     * all particle pairs in the sub-cell, since then it is likely
     * that the box pair has atom pairs within the cut-off.
     * We use the nblist cut-off minus 0.5 times the average x/y diagonal
     * spacing of the sub-cells. Around 40% of the checked pairs are pruned.
     * Using more than 0.5 gains at most 0.5%.
     * If forces are calculated more than twice, the performance gain
     * in the force calculation outweighs the cost of checking.
     * Note that with subcell lists, the atom-pair distance check
     * is only performed when only 1 out of 8 sub-cells in within range,
     * this is because the GPU is much faster than the cpu.
     */
    real bbx, bby;
    real rbb2;

    bbx = 0.5*(gridi->sx + gridj->sx);
    bby = 0.5*(gridi->sy + gridj->sy);
    if (!simple)
    {
        bbx /= GPU_NSUBCELL_X;
        bby /= GPU_NSUBCELL_Y;
    }

    rbb2 = sqr(max(0, rlist - 0.5*sqrt(bbx*bbx + bby*bby)));

#ifndef GMX_DOUBLE
    return rbb2;
#else
    return (float)((1+GMX_FLOAT_EPS)*rbb2);
#endif
}

static int get_ci_block_size(const nbnxn_grid_t *gridi,
                             gmx_bool bDomDec, int nth)
{
    const int ci_block_enum      = 5;
    const int ci_block_denom     = 11;
    const int ci_block_min_atoms = 16;
    int       ci_block;

    /* Here we decide how to distribute the blocks over the threads.
     * We use prime numbers to try to avoid that the grid size becomes
     * a multiple of the number of threads, which would lead to some
     * threads getting "inner" pairs and others getting boundary pairs,
     * which in turns will lead to load imbalance between threads.
     * Set the block size as 5/11/ntask times the average number of cells
     * in a y,z slab. This should ensure a quite uniform distribution
     * of the grid parts of the different thread along all three grid
     * zone boundaries with 3D domain decomposition. At the same time
     * the blocks will not become too small.
     */
    ci_block = (gridi->nc*ci_block_enum)/(ci_block_denom*gridi->ncx*nth);

    /* Ensure the blocks are not too small: avoids cache invalidation */
    if (ci_block*gridi->na_sc < ci_block_min_atoms)
    {
        ci_block = (ci_block_min_atoms + gridi->na_sc - 1)/gridi->na_sc;
    }

    /* Without domain decomposition
     * or with less than 3 blocks per task, divide in nth blocks.
     */
    if (!bDomDec || nth*3*ci_block > gridi->nc)
    {
        ci_block = (gridi->nc + nth - 1)/nth;
    }

    if (ci_block > 1 && (nth - 1)*ci_block >= gridi->nc)
    {
        /* Some threads have no work. Although reducing the block size
         * does not decrease the block count on the first few threads,
         * with GPUs better mixing of "upper" cells that have more empty
         * clusters results in a somewhat lower max load over all threads.
         * Without GPUs the regime of so few atoms per thread is less
         * performance relevant, but with 8-wide SIMD the same reasoning
         * applies, since the pair list uses 4 i-atom "sub-clusters".
         */
        ci_block--;
    }

    return ci_block;
}

/* Generates the part of pair-list nbl assigned to our thread */
static void nbnxn_make_pairlist_part(const nbnxn_search_t nbs,
                                     const nbnxn_grid_t *gridi,
                                     const nbnxn_grid_t *gridj,
                                     nbnxn_search_work_t *work,
                                     const nbnxn_atomdata_t *nbat,
                                     const t_blocka *excl,
                                     real rlist,
                                     int nb_kernel_type,
                                     int ci_block,
                                     gmx_bool bFBufferFlag,
                                     int nsubpair_max,
                                     gmx_bool progBal,
                                     int nsubpair_tot_est,
                                     int th, int nth,
                                     nbnxn_pairlist_t *nbl,
                                     t_nblist *nbl_fep)
{
    int               na_cj_2log;
    matrix            box;
    real              rl2, rl_fep2 = 0;
    float             rbb2;
    int               d;
    int               ci_b, ci, ci_x, ci_y, ci_xy, cj;
    ivec              shp;
    int               tx, ty, tz;
    int               shift;
    gmx_bool          bMakeList;
    real              shx, shy, shz;
    int               conv_i, cell0_i;
    const nbnxn_bb_t *bb_i = NULL;
#ifdef NBNXN_BBXXXX
    const float      *pbb_i = NULL;
#endif
    const float      *bbcz_i, *bbcz_j;
    const int        *flags_i;
    real              bx0, bx1, by0, by1, bz0, bz1;
    real              bz1_frac;
    real              d2cx, d2z, d2z_cx, d2z_cy, d2zx, d2zxy, d2xy;
    int               cxf, cxl, cyf, cyf_x, cyl;
    int               cx, cy;
    int               c0, c1, cs, cf, cl;
    int               ndistc;
    int               ncpcheck;
    int               gridi_flag_shift = 0, gridj_flag_shift = 0;
    gmx_bitmask_t    *gridj_flag       = NULL;
    int               ncj_old_i, ncj_old_j;

    nbs_cycle_start(&work->cc[enbsCCsearch]);

    if (gridj->bSimple != nbl->bSimple)
    {
        gmx_incons("Grid incompatible with pair-list");
    }

    sync_work(nbl);
    nbl->na_sc = gridj->na_sc;
    nbl->na_ci = gridj->na_c;
    nbl->na_cj = nbnxn_kernel_to_cj_size(nb_kernel_type);
    na_cj_2log = get_2log(nbl->na_cj);

    nbl->rlist  = rlist;

    if (bFBufferFlag)
    {
        /* Determine conversion of clusters to flag blocks */
        gridi_flag_shift = 0;
        while ((nbl->na_ci<<gridi_flag_shift) < NBNXN_BUFFERFLAG_SIZE)
        {
            gridi_flag_shift++;
        }
        gridj_flag_shift = 0;
        while ((nbl->na_cj<<gridj_flag_shift) < NBNXN_BUFFERFLAG_SIZE)
        {
            gridj_flag_shift++;
        }

        gridj_flag = work->buffer_flags.flag;
    }

    copy_mat(nbs->box, box);

    rl2 = nbl->rlist*nbl->rlist;

    if (nbs->bFEP && !nbl->bSimple)
    {
        /* Determine an atom-pair list cut-off distance for FEP atom pairs.
         * We should not simply use rlist, since then we would not have
         * the small, effective buffering of the NxN lists.
         * The buffer is on overestimate, but the resulting cost for pairs
         * beyond rlist is neglible compared to the FEP pairs within rlist.
         */
        rl_fep2 = nbl->rlist + effective_buffer_1x1_vs_MxN(gridi, gridj);

        if (debug)
        {
            fprintf(debug, "nbl_fep atom-pair rlist %f\n", rl_fep2);
        }
        rl_fep2 = rl_fep2*rl_fep2;
    }

    rbb2 = boundingbox_only_distance2(gridi, gridj, nbl->rlist, nbl->bSimple);

    if (debug)
    {
        fprintf(debug, "nbl bounding box only distance %f\n", sqrt(rbb2));
    }

    /* Set the shift range */
    for (d = 0; d < DIM; d++)
    {
        /* Check if we need periodicity shifts.
         * Without PBC or with domain decomposition we don't need them.
         */
        if (d >= ePBC2npbcdim(nbs->ePBC) || nbs->dd_dim[d])
        {
            shp[d] = 0;
        }
        else
        {
            if (d == XX &&
                box[XX][XX] - fabs(box[YY][XX]) - fabs(box[ZZ][XX]) < sqrt(rl2))
            {
                shp[d] = 2;
            }
            else
            {
                shp[d] = 1;
            }
        }
    }

    if (nbl->bSimple && !gridi->bSimple)
    {
        conv_i  = gridi->na_sc/gridj->na_sc;
        bb_i    = gridi->bb_simple;
        bbcz_i  = gridi->bbcz_simple;
        flags_i = gridi->flags_simple;
    }
    else
    {
        conv_i  = 1;
#ifdef NBNXN_BBXXXX
        if (gridi->bSimple)
        {
            bb_i  = gridi->bb;
        }
        else
        {
            pbb_i = gridi->pbb;
        }
#else
        /* We use the normal bounding box format for both grid types */
        bb_i  = gridi->bb;
#endif
        bbcz_i  = gridi->bbcz;
        flags_i = gridi->flags;
    }
    cell0_i = gridi->cell0*conv_i;

    bbcz_j = gridj->bbcz;

    if (conv_i != 1)
    {
        /* Blocks of the conversion factor - 1 give a large repeat count
         * combined with a small block size. This should result in good
         * load balancing for both small and large domains.
         */
        ci_block = conv_i - 1;
    }
    if (debug)
    {
        fprintf(debug, "nbl nc_i %d col.av. %.1f ci_block %d\n",
                gridi->nc, gridi->nc/(double)(gridi->ncx*gridi->ncy), ci_block);
    }

    ndistc   = 0;
    ncpcheck = 0;

    /* Initially ci_b and ci to 1 before where we want them to start,
     * as they will both be incremented in next_ci.
     */
    ci_b = -1;
    ci   = th*ci_block - 1;
    ci_x = 0;
    ci_y = 0;
    while (next_ci(gridi, conv_i, nth, ci_block, &ci_x, &ci_y, &ci_b, &ci))
    {
        if (nbl->bSimple && flags_i[ci] == 0)
        {
            continue;
        }

        ncj_old_i = nbl->ncj;

        d2cx = 0;
        if (gridj != gridi && shp[XX] == 0)
        {
            if (nbl->bSimple)
            {
                bx1 = bb_i[ci].upper[BB_X];
            }
            else
            {
                bx1 = gridi->c0[XX] + (ci_x+1)*gridi->sx;
            }
            if (bx1 < gridj->c0[XX])
            {
                d2cx = sqr(gridj->c0[XX] - bx1);

                if (d2cx >= rl2)
                {
                    continue;
                }
            }
        }

        ci_xy = ci_x*gridi->ncy + ci_y;

        /* Loop over shift vectors in three dimensions */
        for (tz = -shp[ZZ]; tz <= shp[ZZ]; tz++)
        {
            shz = tz*box[ZZ][ZZ];

            bz0 = bbcz_i[ci*NNBSBB_D  ] + shz;
            bz1 = bbcz_i[ci*NNBSBB_D+1] + shz;

            if (tz == 0)
            {
                d2z = 0;
            }
            else if (tz < 0)
            {
                d2z = sqr(bz1);
            }
            else
            {
                d2z = sqr(bz0 - box[ZZ][ZZ]);
            }

            d2z_cx = d2z + d2cx;

            if (d2z_cx >= rl2)
            {
                continue;
            }

            bz1_frac =
                bz1/((real)(gridi->cxy_ind[ci_xy+1] - gridi->cxy_ind[ci_xy]));
            if (bz1_frac < 0)
            {
                bz1_frac = 0;
            }
            /* The check with bz1_frac close to or larger than 1 comes later */

            for (ty = -shp[YY]; ty <= shp[YY]; ty++)
            {
                shy = ty*box[YY][YY] + tz*box[ZZ][YY];

                if (nbl->bSimple)
                {
                    by0 = bb_i[ci].lower[BB_Y] + shy;
                    by1 = bb_i[ci].upper[BB_Y] + shy;
                }
                else
                {
                    by0 = gridi->c0[YY] + (ci_y  )*gridi->sy + shy;
                    by1 = gridi->c0[YY] + (ci_y+1)*gridi->sy + shy;
                }

                get_cell_range(by0, by1,
                               gridj->ncy, gridj->c0[YY], gridj->sy, gridj->inv_sy,
                               d2z_cx, rl2,
                               &cyf, &cyl);

                if (cyf > cyl)
                {
                    continue;
                }

                d2z_cy = d2z;
                if (by1 < gridj->c0[YY])
                {
                    d2z_cy += sqr(gridj->c0[YY] - by1);
                }
                else if (by0 > gridj->c1[YY])
                {
                    d2z_cy += sqr(by0 - gridj->c1[YY]);
                }

                for (tx = -shp[XX]; tx <= shp[XX]; tx++)
                {
                    shift = XYZ2IS(tx, ty, tz);

#ifdef NBNXN_SHIFT_BACKWARD
                    if (gridi == gridj && shift > CENTRAL)
                    {
                        continue;
                    }
#endif

                    shx = tx*box[XX][XX] + ty*box[YY][XX] + tz*box[ZZ][XX];

                    if (nbl->bSimple)
                    {
                        bx0 = bb_i[ci].lower[BB_X] + shx;
                        bx1 = bb_i[ci].upper[BB_X] + shx;
                    }
                    else
                    {
                        bx0 = gridi->c0[XX] + (ci_x  )*gridi->sx + shx;
                        bx1 = gridi->c0[XX] + (ci_x+1)*gridi->sx + shx;
                    }

                    get_cell_range(bx0, bx1,
                                   gridj->ncx, gridj->c0[XX], gridj->sx, gridj->inv_sx,
                                   d2z_cy, rl2,
                                   &cxf, &cxl);

                    if (cxf > cxl)
                    {
                        continue;
                    }

                    if (nbl->bSimple)
                    {
                        new_ci_entry(nbl, cell0_i+ci, shift, flags_i[ci]);
                    }
                    else
                    {
                        new_sci_entry(nbl, cell0_i+ci, shift);
                    }

#ifndef NBNXN_SHIFT_BACKWARD
                    if (cxf < ci_x)
#else
                    if (shift == CENTRAL && gridi == gridj &&
                        cxf < ci_x)
#endif
                    {
                        /* Leave the pairs with i > j.
                         * x is the major index, so skip half of it.
                         */
                        cxf = ci_x;
                    }

                    if (nbl->bSimple)
                    {
                        set_icell_bb_simple(bb_i, ci, shx, shy, shz,
                                            nbl->work->bb_ci);
                    }
                    else
                    {
#ifdef NBNXN_BBXXXX
                        set_icell_bbxxxx_supersub(pbb_i, ci, shx, shy, shz,
                                                  nbl->work->pbb_ci);
#else
                        set_icell_bb_supersub(bb_i, ci, shx, shy, shz,
                                              nbl->work->bb_ci);
#endif
                    }

                    nbs->icell_set_x(cell0_i+ci, shx, shy, shz,
                                     gridi->na_c, nbat->xstride, nbat->x,
                                     nbl->work);

                    for (cx = cxf; cx <= cxl; cx++)
                    {
                        d2zx = d2z;
                        if (gridj->c0[XX] + cx*gridj->sx > bx1)
                        {
                            d2zx += sqr(gridj->c0[XX] + cx*gridj->sx - bx1);
                        }
                        else if (gridj->c0[XX] + (cx+1)*gridj->sx < bx0)
                        {
                            d2zx += sqr(gridj->c0[XX] + (cx+1)*gridj->sx - bx0);
                        }

#ifndef NBNXN_SHIFT_BACKWARD
                        if (gridi == gridj &&
                            cx == 0 && cyf < ci_y)
#else
                        if (gridi == gridj &&
                            cx == 0 && shift == CENTRAL && cyf < ci_y)
#endif
                        {
                            /* Leave the pairs with i > j.
                             * Skip half of y when i and j have the same x.
                             */
                            cyf_x = ci_y;
                        }
                        else
                        {
                            cyf_x = cyf;
                        }

                        for (cy = cyf_x; cy <= cyl; cy++)
                        {
                            c0 = gridj->cxy_ind[cx*gridj->ncy+cy];
                            c1 = gridj->cxy_ind[cx*gridj->ncy+cy+1];
#ifdef NBNXN_SHIFT_BACKWARD
                            if (gridi == gridj &&
                                shift == CENTRAL && c0 < ci)
                            {
                                c0 = ci;
                            }
#endif

                            d2zxy = d2zx;
                            if (gridj->c0[YY] + cy*gridj->sy > by1)
                            {
                                d2zxy += sqr(gridj->c0[YY] + cy*gridj->sy - by1);
                            }
                            else if (gridj->c0[YY] + (cy+1)*gridj->sy < by0)
                            {
                                d2zxy += sqr(gridj->c0[YY] + (cy+1)*gridj->sy - by0);
                            }
                            if (c1 > c0 && d2zxy < rl2)
                            {
                                cs = c0 + (int)(bz1_frac*(c1 - c0));
                                if (cs >= c1)
                                {
                                    cs = c1 - 1;
                                }

                                d2xy = d2zxy - d2z;

                                /* Find the lowest cell that can possibly
                                 * be within range.
                                 */
                                cf = cs;
                                while (cf > c0 &&
                                       (bbcz_j[cf*NNBSBB_D+1] >= bz0 ||
                                        d2xy + sqr(bbcz_j[cf*NNBSBB_D+1] - bz0) < rl2))
                                {
                                    cf--;
                                }

                                /* Find the highest cell that can possibly
                                 * be within range.
                                 */
                                cl = cs;
                                while (cl < c1-1 &&
                                       (bbcz_j[cl*NNBSBB_D] <= bz1 ||
                                        d2xy + sqr(bbcz_j[cl*NNBSBB_D] - bz1) < rl2))
                                {
                                    cl++;
                                }

#ifdef NBNXN_REFCODE
                                {
                                    /* Simple reference code, for debugging,
                                     * overrides the more complex code above.
                                     */
                                    int k;
                                    cf = c1;
                                    cl = -1;
                                    for (k = c0; k < c1; k++)
                                    {
                                        if (box_dist2(bx0, bx1, by0, by1, bz0, bz1, bb+k) < rl2 &&
                                            k < cf)
                                        {
                                            cf = k;
                                        }
                                        if (box_dist2(bx0, bx1, by0, by1, bz0, bz1, bb+k) < rl2 &&
                                            k > cl)
                                        {
                                            cl = k;
                                        }
                                    }
                                }
#endif

                                if (gridi == gridj)
                                {
                                    /* We want each atom/cell pair only once,
                                     * only use cj >= ci.
                                     */
#ifndef NBNXN_SHIFT_BACKWARD
                                    cf = max(cf, ci);
#else
                                    if (shift == CENTRAL)
                                    {
                                        cf = max(cf, ci);
                                    }
#endif
                                }

                                if (cf <= cl)
                                {
                                    /* For f buffer flags with simple lists */
                                    ncj_old_j = nbl->ncj;

                                    switch (nb_kernel_type)
                                    {
                                        case nbnxnk4x4_PlainC:
                                            check_subcell_list_space_simple(nbl, cl-cf+1);

                                            make_cluster_list_simple(gridj,
                                                                     nbl, ci, cf, cl,
                                                                     (gridi == gridj && shift == CENTRAL),
                                                                     nbat->x,
                                                                     rl2, rbb2,
                                                                     &ndistc);
                                            break;
#ifdef GMX_NBNXN_SIMD_4XN
                                        case nbnxnk4xN_SIMD_4xN:
                                            check_subcell_list_space_simple(nbl, ci_to_cj(na_cj_2log, cl-cf)+2);
                                            make_cluster_list_simd_4xn(gridj,
                                                                       nbl, ci, cf, cl,
                                                                       (gridi == gridj && shift == CENTRAL),
                                                                       nbat->x,
                                                                       rl2, rbb2,
                                                                       &ndistc);
                                            break;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
                                        case nbnxnk4xN_SIMD_2xNN:
                                            check_subcell_list_space_simple(nbl, ci_to_cj(na_cj_2log, cl-cf)+2);
                                            make_cluster_list_simd_2xnn(gridj,
                                                                        nbl, ci, cf, cl,
                                                                        (gridi == gridj && shift == CENTRAL),
                                                                        nbat->x,
                                                                        rl2, rbb2,
                                                                        &ndistc);
                                            break;
#endif
                                        case nbnxnk8x8x8_PlainC:
                                        case nbnxnk8x8x8_GPU:
                                            check_subcell_list_space_supersub(nbl, cl-cf+1);
                                            for (cj = cf; cj <= cl; cj++)
                                            {
                                                make_cluster_list_supersub(gridi, gridj,
                                                                           nbl, ci, cj,
                                                                           (gridi == gridj && shift == CENTRAL && ci == cj),
                                                                           nbat->xstride, nbat->x,
                                                                           rl2, rbb2,
                                                                           &ndistc);
                                            }
                                            break;
                                    }
                                    ncpcheck += cl - cf + 1;

                                    if (bFBufferFlag && nbl->ncj > ncj_old_j)
                                    {
                                        int cbf, cbl, cb;

                                        cbf = nbl->cj[ncj_old_j].cj >> gridj_flag_shift;
                                        cbl = nbl->cj[nbl->ncj-1].cj >> gridj_flag_shift;
                                        for (cb = cbf; cb <= cbl; cb++)
                                        {
                                            bitmask_init_bit(&gridj_flag[cb], th);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    /* Set the exclusions for this ci list */
                    if (nbl->bSimple)
                    {
                        set_ci_top_excls(nbs,
                                         nbl,
                                         shift == CENTRAL && gridi == gridj,
                                         gridj->na_c_2log,
                                         na_cj_2log,
                                         &(nbl->ci[nbl->nci]),
                                         excl);

                        if (nbs->bFEP)
                        {
                            make_fep_list(nbs, nbat, nbl,
                                          shift == CENTRAL && gridi == gridj,
                                          &(nbl->ci[nbl->nci]),
                                          gridi, gridj, nbl_fep);
                        }
                    }
                    else
                    {
                        set_sci_top_excls(nbs,
                                          nbl,
                                          shift == CENTRAL && gridi == gridj,
                                          gridj->na_c_2log,
                                          &(nbl->sci[nbl->nsci]),
                                          excl);

                        if (nbs->bFEP)
                        {
                            make_fep_list_supersub(nbs, nbat, nbl,
                                                   shift == CENTRAL && gridi == gridj,
                                                   &(nbl->sci[nbl->nsci]),
                                                   shx, shy, shz,
                                                   rl_fep2,
                                                   gridi, gridj, nbl_fep);
                        }
                    }

                    /* Close this ci list */
                    if (nbl->bSimple)
                    {
                        close_ci_entry_simple(nbl);
                    }
                    else
                    {
                        close_ci_entry_supersub(nbl,
                                                nsubpair_max,
                                                progBal, nsubpair_tot_est,
                                                th, nth);
                    }
                }
            }
        }

        if (bFBufferFlag && nbl->ncj > ncj_old_i)
        {
            bitmask_init_bit(&(work->buffer_flags.flag[(gridi->cell0+ci)>>gridi_flag_shift]), th);
        }
    }

    work->ndistc = ndistc;

    nbs_cycle_stop(&work->cc[enbsCCsearch]);

    if (debug)
    {
        fprintf(debug, "number of distance checks %d\n", ndistc);
        fprintf(debug, "ncpcheck %s %d\n", gridi == gridj ? "local" : "non-local",
                ncpcheck);

        if (nbl->bSimple)
        {
            print_nblist_statistics_simple(debug, nbl, nbs, rlist);
        }
        else
        {
            print_nblist_statistics_supersub(debug, nbl, nbs, rlist);
        }

        if (nbs->bFEP)
        {
            fprintf(debug, "nbl FEP list pairs: %d\n", nbl_fep->nrj);
        }
    }
}

static void reduce_buffer_flags(const nbnxn_search_t        nbs,
                                int                         nsrc,
                                const nbnxn_buffer_flags_t *dest)
{
    int            s, b;
    gmx_bitmask_t *flag;

    for (s = 0; s < nsrc; s++)
    {
        flag = nbs->work[s].buffer_flags.flag;

        for (b = 0; b < dest->nflag; b++)
        {
            bitmask_union(&(dest->flag[b]), flag[b]);
        }
    }
}

static void print_reduction_cost(const nbnxn_buffer_flags_t *flags, int nout)
{
    int           nelem, nkeep, ncopy, nred, b, c, out;
    gmx_bitmask_t mask_0;

    nelem = 0;
    nkeep = 0;
    ncopy = 0;
    nred  = 0;
    bitmask_init_bit(&mask_0, 0);
    for (b = 0; b < flags->nflag; b++)
    {
        if (bitmask_is_equal(flags->flag[b], mask_0))
        {
            /* Only flag 0 is set, no copy of reduction required */
            nelem++;
            nkeep++;
        }
        else if (!bitmask_is_zero(flags->flag[b]))
        {
            c = 0;
            for (out = 0; out < nout; out++)
            {
                if (bitmask_is_set(flags->flag[b], out))
                {
                    c++;
                }
            }
            nelem += c;
            if (c == 1)
            {
                ncopy++;
            }
            else
            {
                nred += c;
            }
        }
    }

    fprintf(debug, "nbnxn reduction: #flag %d #list %d elem %4.2f, keep %4.2f copy %4.2f red %4.2f\n",
            flags->nflag, nout,
            nelem/(double)(flags->nflag),
            nkeep/(double)(flags->nflag),
            ncopy/(double)(flags->nflag),
            nred/(double)(flags->nflag));
}

/* Perform a count (linear) sort to sort the smaller lists to the end.
 * This avoids load imbalance on the GPU, as large lists will be
 * scheduled and executed first and the smaller lists later.
 * Load balancing between multi-processors only happens at the end
 * and there smaller lists lead to more effective load balancing.
 * The sorting is done on the cj4 count, not on the actual pair counts.
 * Not only does this make the sort faster, but it also results in
 * better load balancing than using a list sorted on exact load.
 * This function swaps the pointer in the pair list to avoid a copy operation.
 */
static void sort_sci(nbnxn_pairlist_t *nbl)
{
    nbnxn_list_work_t *work;
    int                m, i, s, s0, s1;
    nbnxn_sci_t       *sci_sort;

    if (nbl->ncj4 <= nbl->nsci)
    {
        /* nsci = 0 or all sci have size 1, sorting won't change the order */
        return;
    }

    work = nbl->work;

    /* We will distinguish differences up to double the average */
    m = (2*nbl->ncj4)/nbl->nsci;

    if (m + 1 > work->sort_nalloc)
    {
        work->sort_nalloc = over_alloc_large(m + 1);
        srenew(work->sort, work->sort_nalloc);
    }

    if (work->sci_sort_nalloc != nbl->sci_nalloc)
    {
        work->sci_sort_nalloc = nbl->sci_nalloc;
        nbnxn_realloc_void((void **)&work->sci_sort,
                           0,
                           work->sci_sort_nalloc*sizeof(*work->sci_sort),
                           nbl->alloc, nbl->free);
    }

    /* Count the entries of each size */
    for (i = 0; i <= m; i++)
    {
        work->sort[i] = 0;
    }
    for (s = 0; s < nbl->nsci; s++)
    {
        i = min(m, nbl->sci[s].cj4_ind_end - nbl->sci[s].cj4_ind_start);
        work->sort[i]++;
    }
    /* Calculate the offset for each count */
    s0            = work->sort[m];
    work->sort[m] = 0;
    for (i = m - 1; i >= 0; i--)
    {
        s1            = work->sort[i];
        work->sort[i] = work->sort[i + 1] + s0;
        s0            = s1;
    }

    /* Sort entries directly into place */
    sci_sort = work->sci_sort;
    for (s = 0; s < nbl->nsci; s++)
    {
        i = min(m, nbl->sci[s].cj4_ind_end - nbl->sci[s].cj4_ind_start);
        sci_sort[work->sort[i]++] = nbl->sci[s];
    }

    /* Swap the sci pointers so we use the new, sorted list */
    work->sci_sort = nbl->sci;
    nbl->sci       = sci_sort;
}

/* Make a local or non-local pair-list, depending on iloc */
void nbnxn_make_pairlist(const nbnxn_search_t  nbs,
                         nbnxn_atomdata_t     *nbat,
                         const t_blocka       *excl,
                         real                  rlist,
                         int                   min_ci_balanced,
                         nbnxn_pairlist_set_t *nbl_list,
                         int                   iloc,
                         int                   nb_kernel_type,
                         t_nrnb               *nrnb)
{
    nbnxn_grid_t      *gridi, *gridj;
    gmx_bool           bGPUCPU;
    int                nzi, zi, zj0, zj1, zj;
    int                nsubpair_target, nsubpair_tot_est;
    int                th;
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                ci_block;
    gmx_bool           CombineNBLists;
    gmx_bool           progBal;
    int                np_tot, np_noq, np_hlj, nap;

    /* Check if we are running hybrid GPU + CPU nbnxn mode */
    bGPUCPU = (!nbs->grid[0].bSimple && nbl_list->bSimple);

    nnbl            = nbl_list->nnbl;
    nbl             = nbl_list->nbl;
    CombineNBLists  = nbl_list->bCombined;

    if (debug)
    {
        fprintf(debug, "ns making %d nblists\n", nnbl);
    }

    nbat->bUseBufferFlags = (nbat->nout > 1);
    /* We should re-init the flags before making the first list */
    if (nbat->bUseBufferFlags && (LOCAL_I(iloc) || bGPUCPU))
    {
        init_buffer_flags(&nbat->buffer_flags, nbat->natoms);
    }

    if (nbl_list->bSimple)
    {
        switch (nb_kernel_type)
        {
#ifdef GMX_NBNXN_SIMD_4XN
            case nbnxnk4xN_SIMD_4xN:
                nbs->icell_set_x = icell_set_x_simd_4xn;
                break;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
            case nbnxnk4xN_SIMD_2xNN:
                nbs->icell_set_x = icell_set_x_simd_2xnn;
                break;
#endif
            default:
                nbs->icell_set_x = icell_set_x_simple;
                break;
        }
    }
    else
    {
#ifdef NBNXN_SEARCH_BB_SIMD4
        nbs->icell_set_x = icell_set_x_supersub_simd4;
#else
        nbs->icell_set_x = icell_set_x_supersub;
#endif
    }

    if (LOCAL_I(iloc))
    {
        /* Only zone (grid) 0 vs 0 */
        nzi = 1;
        zj0 = 0;
        zj1 = 1;
    }
    else
    {
        nzi = nbs->zones->nizone;
    }

    if (!nbl_list->bSimple && min_ci_balanced > 0)
    {
        get_nsubpair_target(nbs, iloc, rlist, min_ci_balanced,
                            &nsubpair_target, &nsubpair_tot_est);
    }
    else
    {
        nsubpair_target  = 0;
        nsubpair_tot_est = 0;
    }

    /* Clear all pair-lists */
    for (th = 0; th < nnbl; th++)
    {
        clear_pairlist(nbl[th]);

        if (nbs->bFEP)
        {
            clear_pairlist_fep(nbl_list->nbl_fep[th]);
        }
    }

    for (zi = 0; zi < nzi; zi++)
    {
        gridi = &nbs->grid[zi];

        if (NONLOCAL_I(iloc))
        {
            zj0 = nbs->zones->izone[zi].j0;
            zj1 = nbs->zones->izone[zi].j1;
            if (zi == 0)
            {
                zj0++;
            }
        }
        for (zj = zj0; zj < zj1; zj++)
        {
            gridj = &nbs->grid[zj];

            if (debug)
            {
                fprintf(debug, "ns search grid %d vs %d\n", zi, zj);
            }

            nbs_cycle_start(&nbs->cc[enbsCCsearch]);

            if (nbl[0]->bSimple && !gridi->bSimple)
            {
                /* Hybrid list, determine blocking later */
                ci_block = 0;
            }
            else
            {
                ci_block = get_ci_block_size(gridi, nbs->DomDec, nnbl);
            }

            /* With GPU: generate progressively smaller lists for
             * load balancing for local only or non-local with 2 zones.
             */
            progBal = (LOCAL_I(iloc) || nbs->zones->n <= 2);

#pragma omp parallel for num_threads(nnbl) schedule(static)
            for (th = 0; th < nnbl; th++)
            {
                /* Re-init the thread-local work flag data before making
                 * the first list (not an elegant conditional).
                 */
                if (nbat->bUseBufferFlags && ((zi == 0 && zj == 0) ||
                                              (bGPUCPU && zi == 0 && zj == 1)))
                {
                    init_buffer_flags(&nbs->work[th].buffer_flags, nbat->natoms);
                }

                if (CombineNBLists && th > 0)
                {
                    clear_pairlist(nbl[th]);
                }

                /* Divide the i super cell equally over the nblists */
                nbnxn_make_pairlist_part(nbs, gridi, gridj,
                                         &nbs->work[th], nbat, excl,
                                         rlist,
                                         nb_kernel_type,
                                         ci_block,
                                         nbat->bUseBufferFlags,
                                         nsubpair_target,
                                         progBal, nsubpair_tot_est,
                                         th, nnbl,
                                         nbl[th],
                                         nbl_list->nbl_fep[th]);
            }
            nbs_cycle_stop(&nbs->cc[enbsCCsearch]);

            np_tot = 0;
            np_noq = 0;
            np_hlj = 0;
            for (th = 0; th < nnbl; th++)
            {
                inc_nrnb(nrnb, eNR_NBNXN_DIST2, nbs->work[th].ndistc);

                if (nbl_list->bSimple)
                {
                    np_tot += nbl[th]->ncj;
                    np_noq += nbl[th]->work->ncj_noq;
                    np_hlj += nbl[th]->work->ncj_hlj;
                }
                else
                {
                    /* This count ignores potential subsequent pair pruning */
                    np_tot += nbl[th]->nci_tot;
                }
            }
            nap                   = nbl[0]->na_ci*nbl[0]->na_cj;
            nbl_list->natpair_ljq = (np_tot - np_noq)*nap - np_hlj*nap/2;
            nbl_list->natpair_lj  = np_noq*nap;
            nbl_list->natpair_q   = np_hlj*nap/2;

            if (CombineNBLists && nnbl > 1)
            {
                nbs_cycle_start(&nbs->cc[enbsCCcombine]);

                combine_nblists(nnbl-1, nbl+1, nbl[0]);

                nbs_cycle_stop(&nbs->cc[enbsCCcombine]);
            }
        }
    }

    if (!nbl_list->bSimple)
    {
        /* Sort the entries on size, large ones first */
        if (CombineNBLists || nnbl == 1)
        {
            sort_sci(nbl[0]);
        }
        else
        {
#pragma omp parallel for num_threads(nnbl) schedule(static)
            for (th = 0; th < nnbl; th++)
            {
                sort_sci(nbl[th]);
            }
        }
    }

    if (nbat->bUseBufferFlags)
    {
        reduce_buffer_flags(nbs, nnbl, &nbat->buffer_flags);
    }

    if (nbs->bFEP)
    {
        /* Balance the free-energy lists over all the threads */
        balance_fep_lists(nbs, nbl_list);
    }

    /* Special performance logging stuff (env.var. GMX_NBNXN_CYCLE) */
    if (LOCAL_I(iloc))
    {
        nbs->search_count++;
    }
    if (nbs->print_cycles &&
        (!nbs->DomDec || (nbs->DomDec && !LOCAL_I(iloc))) &&
        nbs->search_count % 100 == 0)
    {
        nbs_cycle_print(stderr, nbs);
    }

    if (debug && (CombineNBLists && nnbl > 1))
    {
        if (nbl[0]->bSimple)
        {
            print_nblist_statistics_simple(debug, nbl[0], nbs, rlist);
        }
        else
        {
            print_nblist_statistics_supersub(debug, nbl[0], nbs, rlist);
        }
    }

    if (debug)
    {
        if (gmx_debug_at)
        {
            if (nbl[0]->bSimple)
            {
                print_nblist_ci_cj(debug, nbl[0]);
            }
            else
            {
                print_nblist_sci_cj(debug, nbl[0]);
            }
        }

        if (nbat->bUseBufferFlags)
        {
            print_reduction_cost(&nbat->buffer_flags, nnbl);
        }
    }
}
