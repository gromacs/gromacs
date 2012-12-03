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
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "pbc.h"
#include "nbnxn_consts.h"
#include "nbnxn_internal.h"
#include "nbnxn_atomdata.h"
#include "nbnxn_search.h"
#include "gmx_cyclecounter.h"
#include "gmxfio.h"
#include "gmx_omp_nthreads.h"
#include "nrnb.h"


/* Pair search box lower and upper corner in x,y,z.
 * Store this in 4 iso 3 reals, which is useful with SSE.
 * To avoid complicating the code we also use 4 without SSE.
 */
#define NNBSBB_C         4
#define NNBSBB_B         (2*NNBSBB_C)
/* Pair search box lower and upper bound in z only. */
#define NNBSBB_D         2
/* Pair search box lower and upper corner x,y,z indices */
#define BBL_X  0
#define BBL_Y  1
#define BBL_Z  2
#define BBU_X  4
#define BBU_Y  5
#define BBU_Z  6


#ifdef NBNXN_SEARCH_SSE

#ifndef GMX_DOUBLE
#define NBNXN_SEARCH_SSE_SINGLE
#include "gmx_x86_simd_single.h"
#else
#include "gmx_x86_simd_double.h"
#endif

#if defined NBNXN_SEARCH_SSE_SINGLE && GPU_NSUBCELL == 8
#define NBNXN_8BB_SSE
#endif

/* The width of SSE/AVX128 with single precision for bounding boxes with GPU.
 * Here AVX-256 turns out to be slightly slower than AVX-128.
 */
#define STRIDE_8BB        4
#define STRIDE_8BB_2LOG   2


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
#ifndef GMX_DOUBLE
/* 128 bits can hold 4 floats */
#define CI_TO_CJ_S128(ci)  CI_TO_CJ_J4(ci)
#define X_IND_CI_S128(ci)  X_IND_CI_J4(ci)
#define X_IND_CJ_S128(cj)  X_IND_CJ_J4(cj)
/* 256 bits can hold 8 floats */
#define CI_TO_CJ_S256(ci)  CI_TO_CJ_J8(ci)
#define X_IND_CI_S256(ci)  X_IND_CI_J8(ci)
#define X_IND_CJ_S256(cj)  X_IND_CJ_J8(cj)
#else
/* 128 bits can hold 2 doubles */
#define CI_TO_CJ_S128(ci)  CI_TO_CJ_J2(ci)
#define X_IND_CI_S128(ci)  X_IND_CI_J2(ci)
#define X_IND_CJ_S128(cj)  X_IND_CJ_J2(cj)
/* 256 bits can hold 4 doubles */
#define CI_TO_CJ_S256(ci)  CI_TO_CJ_J4(ci)
#define X_IND_CI_S256(ci)  X_IND_CI_J4(ci)
#define X_IND_CJ_S256(cj)  X_IND_CJ_J4(cj)
#endif

#endif /* NBNXN_SEARCH_SSE */


/* Interaction masks for 4xN atom interactions.
 * Bit i*CJ_SIZE + j tells if atom i and j interact.
 */
/* All interaction mask is the same for all kernels */
#define NBNXN_INT_MASK_ALL        0xffffffff
/* 4x4 kernel diagonal mask */
#define NBNXN_INT_MASK_DIAG       0x08ce
/* 4x2 kernel diagonal masks */
#define NBNXN_INT_MASK_DIAG_J2_0  0x0002
#define NBNXN_INT_MASK_DIAG_J2_1  0x002F
/* 4x8 kernel diagonal masks */
#define NBNXN_INT_MASK_DIAG_J8_0  0xf0f8fcfe
#define NBNXN_INT_MASK_DIAG_J8_1  0x0080c0e0


#ifdef NBNXN_SEARCH_SSE
/* Store bounding boxes corners as quadruplets: xxxxyyyyzzzz */
#define NBNXN_BBXXXX
/* Size of bounding box corners quadruplet */
#define NNBSBB_XXXX      (NNBSBB_D*DIM*STRIDE_8BB)
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

    for(i=0; i<enbsCCnr; i++)
    {
        cc[i].count = 0;
        cc[i].c     = 0;
    }
}

static double Mcyc_av(const nbnxn_cycle_t *cc)
{
    return (double)cc->c*1e-6/cc->count;
}

static void nbs_cycle_print(FILE *fp,const nbnxn_search_t nbs)
{
    int n;
    int t;

    fprintf(fp,"\n");
    fprintf(fp,"ns %4d grid %4.1f search %4.1f red.f %5.3f",
            nbs->cc[enbsCCgrid].count,
            Mcyc_av(&nbs->cc[enbsCCgrid]),
            Mcyc_av(&nbs->cc[enbsCCsearch]),
            Mcyc_av(&nbs->cc[enbsCCreducef]));

    if (nbs->nthread_max > 1)
    {
        if (nbs->cc[enbsCCcombine].count > 0)
        {
            fprintf(fp," comb %5.2f",
                    Mcyc_av(&nbs->cc[enbsCCcombine]));
        }
        fprintf(fp," s. th");
        for(t=0; t<nbs->nthread_max; t++)
        {
            fprintf(fp," %4.1f",
                    Mcyc_av(&nbs->work[t].cc[enbsCCsearch]));
        }
    }
    fprintf(fp,"\n");
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
        gmx_fatal(FARGS,"nbnxn na_c (%d) is not a power of 2",n);
    }

    return log2;
}

static int nbnxn_kernel_to_ci_size(int nb_kernel_type)
{
    switch (nb_kernel_type)
    {
    case nbk4x4_PlainC:
    case nbk4xN_X86_SIMD128:
    case nbk4xN_X86_SIMD256:
        return NBNXN_CPU_CLUSTER_I_SIZE;
    case nbk8x8x8_CUDA:
    case nbk8x8x8_PlainC:
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
    switch (nb_kernel_type)
    {
    case nbk4x4_PlainC:
        return NBNXN_CPU_CLUSTER_I_SIZE;
    case nbk4xN_X86_SIMD128:
        /* Number of reals that fit in SIMD (128 bits = 16 bytes) */
        return 16/sizeof(real);
    case nbk4xN_X86_SIMD256:
        /* Number of reals that fit in SIMD (256 bits = 32 bytes) */
        return 32/sizeof(real);
    case nbk8x8x8_CUDA:
    case nbk8x8x8_PlainC:
        return nbnxn_kernel_to_ci_size(nb_kernel_type);
    default:
        gmx_incons("unknown kernel type");
    }

    return 0;
}

static int ci_to_cj(int na_cj_2log,int ci)
{
    switch (na_cj_2log)
    {
    case 2: return  ci;     break;
    case 1: return (ci<<1); break;
    case 3: return (ci>>1); break;
    }

    return 0;
}

gmx_bool nbnxn_kernel_pairlist_simple(int nb_kernel_type)
{
    if (nb_kernel_type == nbkNotSet)
    {
        gmx_fatal(FARGS, "Non-bonded kernel type not set for Verlet-style pair-list.");
    }

    switch (nb_kernel_type)
    {
    case nbk8x8x8_CUDA:
    case nbk8x8x8_PlainC:
        return FALSE;

    case nbk4x4_PlainC:
    case nbk4xN_X86_SIMD128:
    case nbk4xN_X86_SIMD256:
        return TRUE;

    default:
        gmx_incons("Invalid nonbonded kernel type passed!");
        return FALSE;
    }
}

void nbnxn_init_search(nbnxn_search_t * nbs_ptr,
                       ivec *n_dd_cells,
                       gmx_domdec_zones_t *zones,
                       int nthread_max)
{
    nbnxn_search_t nbs;
    int d,g,t;

    snew(nbs,1);
    *nbs_ptr = nbs;

    nbs->DomDec = (n_dd_cells != NULL);

    clear_ivec(nbs->dd_dim);
    nbs->ngrid = 1;
    if (nbs->DomDec)
    {
        nbs->zones = zones;

        for(d=0; d<DIM; d++)
        {
            if ((*n_dd_cells)[d] > 1)
            {
                nbs->dd_dim[d] = 1;
                /* Each grid matches a DD zone */
                nbs->ngrid *= 2;
            }
        }
    }

    snew(nbs->grid,nbs->ngrid);
    for(g=0; g<nbs->ngrid; g++)
    {
        nbnxn_grid_init(&nbs->grid[g]);
    }
    nbs->cell        = NULL;
    nbs->cell_nalloc = 0;
    nbs->a           = NULL;
    nbs->a_nalloc    = 0;

    nbs->nthread_max = nthread_max;

    /* Initialize the work data structures for each thread */
    snew(nbs->work,nbs->nthread_max);
    for(t=0; t<nbs->nthread_max; t++)
    {
        nbs->work[t].cxy_na           = NULL;
        nbs->work[t].cxy_na_nalloc    = 0;
        nbs->work[t].sort_work        = NULL;
        nbs->work[t].sort_work_nalloc = 0;
    }

    /* Initialize detailed nbsearch cycle counting */
    nbs->print_cycles = (getenv("GMX_NBNXN_CYCLE") != 0);
    nbs->search_count = 0;
    nbs_cycle_clear(nbs->cc);
    for(t=0; t<nbs->nthread_max; t++)
    {
        nbs_cycle_clear(nbs->work[t].cc);
    }
}

static real grid_atom_density(int n,rvec corner0,rvec corner1)
{
    rvec size;

    rvec_sub(corner1,corner0,size);

    return n/(size[XX]*size[YY]*size[ZZ]);
}

static int set_grid_size_xy(const nbnxn_search_t nbs,
                            nbnxn_grid_t *grid,
                            int n,rvec corner0,rvec corner1,
                            real atom_density,
                            int XFormat)
{
    rvec size;
    int  na_c;
    real adens,tlen,tlen_x,tlen_y,nc_max;
    int  t;

    rvec_sub(corner1,corner0,size);

    if (n > grid->na_sc)
    {
        /* target cell length */
        if (grid->bSimple)
        {
            /* To minimize the zero interactions, we should make
             * the largest of the i/j cell cubic.
             */
            na_c = max(grid->na_c,grid->na_cj);

            /* Approximately cubic cells */
            tlen   = pow(na_c/atom_density,1.0/3.0);
            tlen_x = tlen;
            tlen_y = tlen;
        }
        else
        {
            /* Approximately cubic sub cells */
            tlen   = pow(grid->na_c/atom_density,1.0/3.0);
            tlen_x = tlen*GPU_NSUBCELL_X;
            tlen_y = tlen*GPU_NSUBCELL_Y;
        }
        /* We round ncx and ncy down, because we get less cell pairs
         * in the nbsist when the fixed cell dimensions (x,y) are
         * larger than the variable one (z) than the other way around.
         */
        grid->ncx = max(1,(int)(size[XX]/tlen_x));
        grid->ncy = max(1,(int)(size[YY]/tlen_y));
    }
    else
    {
        grid->ncx = 1;
        grid->ncy = 1;
    }

    /* We need one additional cell entry for particles moved by DD */
    if (grid->ncx*grid->ncy+1 > grid->cxy_nalloc)
    {
        grid->cxy_nalloc = over_alloc_large(grid->ncx*grid->ncy+1);
        srenew(grid->cxy_na,grid->cxy_nalloc);
        srenew(grid->cxy_ind,grid->cxy_nalloc+1);
    }
    for(t=0; t<nbs->nthread_max; t++)
    {
        if (grid->ncx*grid->ncy+1 > nbs->work[t].cxy_na_nalloc)
        {
            nbs->work[t].cxy_na_nalloc = over_alloc_large(grid->ncx*grid->ncy+1);
            srenew(nbs->work[t].cxy_na,nbs->work[t].cxy_na_nalloc);
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
        int bb_nalloc;

        grid->nc_nalloc = over_alloc_large(nc_max);
        srenew(grid->nsubc,grid->nc_nalloc);
        srenew(grid->bbcz,grid->nc_nalloc*NNBSBB_D);
#ifdef NBNXN_8BB_SSE
        bb_nalloc = grid->nc_nalloc*GPU_NSUBCELL/STRIDE_8BB*NNBSBB_XXXX;
#else
        bb_nalloc = grid->nc_nalloc*GPU_NSUBCELL*NNBSBB_B;
#endif
        sfree_aligned(grid->bb);
        /* This snew also zeros the contents, this avoid possible
         * floating exceptions in SSE with the unused bb elements.
         */
        snew_aligned(grid->bb,bb_nalloc,16);

        if (grid->bSimple)
        {
            if (grid->na_cj == grid->na_c)
            {
                grid->bbj = grid->bb;
            }
            else
            {
                sfree_aligned(grid->bbj);
                snew_aligned(grid->bbj,bb_nalloc*grid->na_c/grid->na_cj,16);
            }
        }

        srenew(grid->flags,grid->nc_nalloc);
    }

    copy_rvec(corner0,grid->c0);
    copy_rvec(corner1,grid->c1);
    grid->sx = size[XX]/grid->ncx;
    grid->sy = size[YY]/grid->ncy;
    grid->inv_sx = 1/grid->sx;
    grid->inv_sy = 1/grid->sy;

    return nc_max;
}

#define SORT_GRID_OVERSIZE 2
#define SGSF (SORT_GRID_OVERSIZE + 1)

static void sort_atoms(int dim,gmx_bool Backwards,
                       int *a,int n,rvec *x,
                       real h0,real invh,int nsort,int *sort)
{
    int i,c;
    int zi,zim;
    int cp,tmp;

    if (n <= 1)
    {
        /* Nothing to do */
        return;
    }

    /* For small oversize factors clearing the whole area is fastest.
     * For large oversize we should clear the used elements after use.
     */
    for(i=0; i<nsort; i++)
    {
        sort[i] = -1;
    }
    /* Sort the particles using a simple index sort */
    for(i=0; i<n; i++)
    {
        /* The cast takes care of float-point rounding effects below zero.
         * This code assumes particles are less than 1/SORT_GRID_OVERSIZE
         * times the box height out of the box.
         */
        zi = (int)((x[a[i]][dim] - h0)*invh);

#ifdef DEBUG_NBNXN_GRIDDING
        if (zi < 0 || zi >= nsort)
        {
            gmx_fatal(FARGS,"(int)((x[%d][%c]=%f - %f)*%f) = %d, not in 0 - %d\n",
                      a[i],'x'+dim,x[a[i]][dim],h0,invh,zi,nsort);
        }
#endif

        /* Ideally this particle should go in sort cell zi,
         * but that might already be in use,
         * in that case find the first empty cell higher up
         */
        if (sort[zi] < 0)
        {
            sort[zi] = a[i];
        }
        else
        {
            /* We have multiple atoms in the same sorting slot.
             * Sort on real z for minimal bounding box size.
             * There is an extra check for identical z to ensure
             * well-defined output order, independent of input order
             * to ensure binary reproducibility after restarts.
             */
            while(sort[zi] >= 0 && ( x[a[i]][dim] >  x[sort[zi]][dim] ||
                                    (x[a[i]][dim] == x[sort[zi]][dim] &&
                                     a[i] > sort[zi])))
            {
                zi++;
            }

            if (sort[zi] >= 0)
            {
                /* Shift all elements by one slot until we find an empty slot */
                cp = sort[zi];
                zim = zi + 1;
                while (sort[zim] >= 0)
                {
                    tmp = sort[zim];
                    sort[zim] = cp;
                    cp  = tmp;
                    zim++;
                }
                sort[zim] = cp;
            }
            sort[zi] = a[i];
        }
    }

    c = 0;
    if (!Backwards)
    {
        for(zi=0; zi<nsort; zi++)
        {
            if (sort[zi] >= 0)
            {
                a[c++] = sort[zi];
            }
        }
    }
    else
    {
        for(zi=nsort-1; zi>=0; zi--)
        {
            if (sort[zi] >= 0)
            {
                a[c++] = sort[zi];
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
static void calc_bounding_box(int na,int stride,const real *x,float *bb)
{
    int  i,j;
    real xl,xh,yl,yh,zl,zh;

    i = 0;
    xl = x[i+XX];
    xh = x[i+XX];
    yl = x[i+YY];
    yh = x[i+YY];
    zl = x[i+ZZ];
    zh = x[i+ZZ];
    i += stride;
    for(j=1; j<na; j++)
    {
        xl = min(xl,x[i+XX]);
        xh = max(xh,x[i+XX]);
        yl = min(yl,x[i+YY]);
        yh = max(yh,x[i+YY]);
        zl = min(zl,x[i+ZZ]);
        zh = max(zh,x[i+ZZ]);
        i += stride;
    }
    /* Note: possible double to float conversion here */
    bb[BBL_X] = R2F_D(xl);
    bb[BBL_Y] = R2F_D(yl);
    bb[BBL_Z] = R2F_D(zl);
    bb[BBU_X] = R2F_U(xh);
    bb[BBU_Y] = R2F_U(yh);
    bb[BBU_Z] = R2F_U(zh);
}

/* Packed coordinates, bb order xyz0 */
static void calc_bounding_box_x_x4(int na,const real *x,float *bb)
{
    int  j;
    real xl,xh,yl,yh,zl,zh;

    xl = x[XX*PACK_X4];
    xh = x[XX*PACK_X4];
    yl = x[YY*PACK_X4];
    yh = x[YY*PACK_X4];
    zl = x[ZZ*PACK_X4];
    zh = x[ZZ*PACK_X4];
    for(j=1; j<na; j++)
    {
        xl = min(xl,x[j+XX*PACK_X4]);
        xh = max(xh,x[j+XX*PACK_X4]);
        yl = min(yl,x[j+YY*PACK_X4]);
        yh = max(yh,x[j+YY*PACK_X4]);
        zl = min(zl,x[j+ZZ*PACK_X4]);
        zh = max(zh,x[j+ZZ*PACK_X4]);
    }
    /* Note: possible double to float conversion here */
    bb[BBL_X] = R2F_D(xl);
    bb[BBL_Y] = R2F_D(yl);
    bb[BBL_Z] = R2F_D(zl);
    bb[BBU_X] = R2F_U(xh);
    bb[BBU_Y] = R2F_U(yh);
    bb[BBU_Z] = R2F_U(zh);
}

/* Packed coordinates, bb order xyz0 */
static void calc_bounding_box_x_x8(int na,const real *x,float *bb)
{
    int  j;
    real xl,xh,yl,yh,zl,zh;

    xl = x[XX*PACK_X8];
    xh = x[XX*PACK_X8];
    yl = x[YY*PACK_X8];
    yh = x[YY*PACK_X8];
    zl = x[ZZ*PACK_X8];
    zh = x[ZZ*PACK_X8];
    for(j=1; j<na; j++)
    {
        xl = min(xl,x[j+XX*PACK_X8]);
        xh = max(xh,x[j+XX*PACK_X8]);
        yl = min(yl,x[j+YY*PACK_X8]);
        yh = max(yh,x[j+YY*PACK_X8]);
        zl = min(zl,x[j+ZZ*PACK_X8]);
        zh = max(zh,x[j+ZZ*PACK_X8]);
    }
    /* Note: possible double to float conversion here */
    bb[BBL_X] = R2F_D(xl);
    bb[BBL_Y] = R2F_D(yl);
    bb[BBL_Z] = R2F_D(zl);
    bb[BBU_X] = R2F_U(xh);
    bb[BBU_Y] = R2F_U(yh);
    bb[BBU_Z] = R2F_U(zh);
}

#ifdef NBNXN_SEARCH_SSE

/* Packed coordinates, bb order xyz0 */
static void calc_bounding_box_x_x4_halves(int na,const real *x,
                                          float *bb,float *bbj)
{
    calc_bounding_box_x_x4(min(na,2),x,bbj);

    if (na > 2)
    {
        calc_bounding_box_x_x4(min(na-2,2),x+(PACK_X4>>1),bbj+NNBSBB_B);
    }
    else
    {
        /* Set the "empty" bounding box to the same as the first one,
         * so we don't need to treat special cases in the rest of the code.
         */
        _mm_store_ps(bbj+NNBSBB_B         ,_mm_load_ps(bbj));
        _mm_store_ps(bbj+NNBSBB_B+NNBSBB_C,_mm_load_ps(bbj+NNBSBB_C));
    }

    _mm_store_ps(bb         ,_mm_min_ps(_mm_load_ps(bbj),
                                        _mm_load_ps(bbj+NNBSBB_B)));
    _mm_store_ps(bb+NNBSBB_C,_mm_max_ps(_mm_load_ps(bbj+NNBSBB_C),
                                        _mm_load_ps(bbj+NNBSBB_B+NNBSBB_C)));
}

/* Coordinate order xyz, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx(int na,int stride,const real *x,float *bb)
{
    int  i,j;
    real xl,xh,yl,yh,zl,zh;

    i = 0;
    xl = x[i+XX];
    xh = x[i+XX];
    yl = x[i+YY];
    yh = x[i+YY];
    zl = x[i+ZZ];
    zh = x[i+ZZ];
    i += stride;
    for(j=1; j<na; j++)
    {
        xl = min(xl,x[i+XX]);
        xh = max(xh,x[i+XX]);
        yl = min(yl,x[i+YY]);
        yh = max(yh,x[i+YY]);
        zl = min(zl,x[i+ZZ]);
        zh = max(zh,x[i+ZZ]);
        i += stride;
    }
    /* Note: possible double to float conversion here */
    bb[0*STRIDE_8BB] = R2F_D(xl);
    bb[1*STRIDE_8BB] = R2F_D(yl);
    bb[2*STRIDE_8BB] = R2F_D(zl);
    bb[3*STRIDE_8BB] = R2F_U(xh);
    bb[4*STRIDE_8BB] = R2F_U(yh);
    bb[5*STRIDE_8BB] = R2F_U(zh);
}

#endif /* NBNXN_SEARCH_SSE */

#ifdef NBNXN_SEARCH_SSE_SINGLE

/* Coordinate order xyz?, bb order xyz0 */
static void calc_bounding_box_sse(int na,const float *x,float *bb)
{
    __m128 bb_0_SSE,bb_1_SSE;
    __m128 x_SSE;

    int  i;

    bb_0_SSE = _mm_load_ps(x);
    bb_1_SSE = bb_0_SSE;

    for(i=1; i<na; i++)
    {
        x_SSE    = _mm_load_ps(x+i*NNBSBB_C);
        bb_0_SSE = _mm_min_ps(bb_0_SSE,x_SSE);
        bb_1_SSE = _mm_max_ps(bb_1_SSE,x_SSE);
    }

    _mm_store_ps(bb  ,bb_0_SSE);
    _mm_store_ps(bb+4,bb_1_SSE);
}

/* Coordinate order xyz?, bb order xxxxyyyyzzzz */
static void calc_bounding_box_xxxx_sse(int na,const float *x,
                                       float *bb_work,
                                       real *bb)
{
    calc_bounding_box_sse(na,x,bb_work);

    bb[0*STRIDE_8BB] = bb_work[BBL_X];
    bb[1*STRIDE_8BB] = bb_work[BBL_Y];
    bb[2*STRIDE_8BB] = bb_work[BBL_Z];
    bb[3*STRIDE_8BB] = bb_work[BBU_X];
    bb[4*STRIDE_8BB] = bb_work[BBU_Y];
    bb[5*STRIDE_8BB] = bb_work[BBU_Z];
}

#endif /* NBNXN_SEARCH_SSE_SINGLE */

#ifdef NBNXN_SEARCH_SSE

/* Combines pairs of consecutive bounding boxes */
static void combine_bounding_box_pairs(nbnxn_grid_t *grid,const float *bb)
{
    int    i,j,sc2,nc2,c2;
    __m128 min_SSE,max_SSE;

    for(i=0; i<grid->ncx*grid->ncy; i++)
    {
        /* Starting bb in a column is expected to be 2-aligned */
        sc2 = grid->cxy_ind[i]>>1;
        /* For odd numbers skip the last bb here */
        nc2 = (grid->cxy_na[i]+3)>>(2+1);
        for(c2=sc2; c2<sc2+nc2; c2++)
        {
            min_SSE = _mm_min_ps(_mm_load_ps(bb+(c2*4+0)*NNBSBB_C),
                                 _mm_load_ps(bb+(c2*4+2)*NNBSBB_C));
            max_SSE = _mm_max_ps(_mm_load_ps(bb+(c2*4+1)*NNBSBB_C),
                                 _mm_load_ps(bb+(c2*4+3)*NNBSBB_C));
            _mm_store_ps(grid->bbj+(c2*2+0)*NNBSBB_C,min_SSE);
            _mm_store_ps(grid->bbj+(c2*2+1)*NNBSBB_C,max_SSE);
        }
        if (((grid->cxy_na[i]+3)>>2) & 1)
        {
            /* Copy the last bb for odd bb count in this column */
            for(j=0; j<NNBSBB_C; j++)
            {
                grid->bbj[(c2*2+0)*NNBSBB_C+j] = bb[(c2*4+0)*NNBSBB_C+j];
                grid->bbj[(c2*2+1)*NNBSBB_C+j] = bb[(c2*4+1)*NNBSBB_C+j];
            }
        }
    }
}

#endif


/* Prints the average bb size, used for debug output */
static void print_bbsizes_simple(FILE *fp,
                                 const nbnxn_search_t nbs,
                                 const nbnxn_grid_t *grid)
{
    int  c,d;
    dvec ba;

    clear_dvec(ba);
    for(c=0; c<grid->nc; c++)
    {
        for(d=0; d<DIM; d++)
        {
            ba[d] += grid->bb[c*NNBSBB_B+NNBSBB_C+d] - grid->bb[c*NNBSBB_B+d];
        }
    }
    dsvmul(1.0/grid->nc,ba,ba);

    fprintf(fp,"ns bb: %4.2f %4.2f %4.2f  %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            nbs->box[XX][XX]/grid->ncx,
            nbs->box[YY][YY]/grid->ncy,
            nbs->box[ZZ][ZZ]*grid->ncx*grid->ncy/grid->nc,
            ba[XX],ba[YY],ba[ZZ],
            ba[XX]*grid->ncx/nbs->box[XX][XX],
            ba[YY]*grid->ncy/nbs->box[YY][YY],
            ba[ZZ]*grid->nc/(grid->ncx*grid->ncy*nbs->box[ZZ][ZZ]));
}

/* Prints the average bb size, used for debug output */
static void print_bbsizes_supersub(FILE *fp,
                                   const nbnxn_search_t nbs,
                                   const nbnxn_grid_t *grid)
{
    int  ns,c,s;
    dvec ba;

    clear_dvec(ba);
    ns = 0;
    for(c=0; c<grid->nc; c++)
    {
#ifdef NBNXN_BBXXXX
        for(s=0; s<grid->nsubc[c]; s+=STRIDE_8BB)
        {
            int cs_w,i,d;

            cs_w = (c*GPU_NSUBCELL + s)/STRIDE_8BB;
            for(i=0; i<STRIDE_8BB; i++)
            {
                for(d=0; d<DIM; d++)
                {
                    ba[d] +=
                        grid->bb[cs_w*NNBSBB_XXXX+(DIM+d)*STRIDE_8BB+i] -
                        grid->bb[cs_w*NNBSBB_XXXX+     d *STRIDE_8BB+i];
                }
            }
        }
#else
        for(s=0; s<grid->nsubc[c]; s++)
        {
            int cs,d;

            cs = c*GPU_NSUBCELL + s;
            for(d=0; d<DIM; d++)
            {
                ba[d] +=
                    grid->bb[cs*NNBSBB_B+NNBSBB_C+d] -
                    grid->bb[cs*NNBSBB_B         +d];
            }
        }
#endif
        ns += grid->nsubc[c];
    }
    dsvmul(1.0/ns,ba,ba);

    fprintf(fp,"ns bb: %4.2f %4.2f %4.2f  %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            nbs->box[XX][XX]/(grid->ncx*GPU_NSUBCELL_X),
            nbs->box[YY][YY]/(grid->ncy*GPU_NSUBCELL_Y),
            nbs->box[ZZ][ZZ]*grid->ncx*grid->ncy/(grid->nc*GPU_NSUBCELL_Z),
            ba[XX],ba[YY],ba[ZZ],
            ba[XX]*grid->ncx*GPU_NSUBCELL_X/nbs->box[XX][XX],
            ba[YY]*grid->ncy*GPU_NSUBCELL_Y/nbs->box[YY][YY],
            ba[ZZ]*grid->nc*GPU_NSUBCELL_Z/(grid->ncx*grid->ncy*nbs->box[ZZ][ZZ]));
}

/* Potentially sorts atoms on LJ coefficients !=0 and ==0.
 * Also sets interaction flags.
 */
void sort_on_lj(nbnxn_atomdata_t *nbat,int na_c,
                int a0,int a1,const int *atinfo,
                int *order,
                int *flags)
{
    int subc,s,a,n1,n2,a_lj_max,i,j;
    int sort1[NBNXN_NA_SC_MAX/GPU_NSUBCELL];
    int sort2[NBNXN_NA_SC_MAX/GPU_NSUBCELL];
    gmx_bool haveQ;

    *flags = 0;

    subc = 0;
    for(s=a0; s<a1; s+=na_c)
    {
        /* Make lists for this (sub-)cell on atoms with and without LJ */
        n1 = 0;
        n2 = 0;
        haveQ = FALSE;
        a_lj_max = -1;
        for(a=s; a<min(s+na_c,a1); a++)
        {
            haveQ = haveQ || GET_CGINFO_HAS_Q(atinfo[order[a]]);

            if (GET_CGINFO_HAS_VDW(atinfo[order[a]]))
            {
                sort1[n1++] = order[a];
                a_lj_max = a;
            }
            else
            {
                sort2[n2++] = order[a];
            }
        }

        /* If we don't have atom with LJ, there's nothing to sort */
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
                    for(i=0; i<n1; i++)
                    {
                        order[a0+i] = sort1[i];
                    }
                    for(j=0; j<n2; j++)
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
               int a0,int a1,
               const int *atinfo,
               rvec *x,
               int sx,int sy, int sz,
               float *bb_work)
{
    int    na,a;
    size_t offset;
    float  *bb_ptr;

    na = a1 - a0;

    if (grid->bSimple)
    {
        sort_on_lj(nbat,grid->na_c,a0,a1,atinfo,nbs->a,
                   grid->flags+(a0>>grid->na_c_2log)-grid->cell0);
    }

    /* Now we have sorted the atoms, set the cell indices */
    for(a=a0; a<a1; a++)
    {
        nbs->cell[nbs->a[a]] = a;
    }

    copy_rvec_to_nbat_real(nbs->a+a0,a1-a0,grid->na_c,x,
                           nbat->XFormat,nbat->x,a0,
                           sx,sy,sz);

    if (nbat->XFormat == nbatX4)
    {
        /* Store the bounding boxes as xyz.xyz. */
        offset = ((a0 - grid->cell0*grid->na_sc)>>grid->na_c_2log)*NNBSBB_B;
        bb_ptr = grid->bb + offset;

#if defined GMX_DOUBLE && defined NBNXN_SEARCH_SSE
        if (2*grid->na_cj == grid->na_c)
        {
            calc_bounding_box_x_x4_halves(na,nbat->x+X4_IND_A(a0),bb_ptr,
                                          grid->bbj+offset*2);
        }
        else
#endif
        {
            calc_bounding_box_x_x4(na,nbat->x+X4_IND_A(a0),bb_ptr);
        }
    }
    else if (nbat->XFormat == nbatX8)
    {
        /* Store the bounding boxes as xyz.xyz. */
        offset = ((a0 - grid->cell0*grid->na_sc)>>grid->na_c_2log)*NNBSBB_B;
        bb_ptr = grid->bb + offset;

        calc_bounding_box_x_x8(na,nbat->x+X8_IND_A(a0),bb_ptr);
    }
#ifdef NBNXN_BBXXXX
    else if (!grid->bSimple)
    {
        /* Store the bounding boxes in a format convenient
         * for SSE calculations: xxxxyyyyzzzz...
                             */
        bb_ptr =
            grid->bb +
            ((a0-grid->cell0*grid->na_sc)>>(grid->na_c_2log+STRIDE_8BB_2LOG))*NNBSBB_XXXX +
            (((a0-grid->cell0*grid->na_sc)>>grid->na_c_2log) & (STRIDE_8BB-1));

#ifdef NBNXN_SEARCH_SSE_SINGLE
        if (nbat->XFormat == nbatXYZQ)
        {
            calc_bounding_box_xxxx_sse(na,nbat->x+a0*nbat->xstride,
                                       bb_work,bb_ptr);
        }
        else
#endif
        {
            calc_bounding_box_xxxx(na,nbat->xstride,nbat->x+a0*nbat->xstride,
                                   bb_ptr);
        }
        if (gmx_debug_at)
        {
            fprintf(debug,"%2d %2d %2d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    sx,sy,sz,
                    bb_ptr[0*STRIDE_8BB],bb_ptr[3*STRIDE_8BB],
                    bb_ptr[1*STRIDE_8BB],bb_ptr[4*STRIDE_8BB],
                    bb_ptr[2*STRIDE_8BB],bb_ptr[5*STRIDE_8BB]);
        }
    }
#endif
    else
    {
        /* Store the bounding boxes as xyz.xyz. */
        bb_ptr = grid->bb+((a0-grid->cell0*grid->na_sc)>>grid->na_c_2log)*NNBSBB_B;

        calc_bounding_box(na,nbat->xstride,nbat->x+a0*nbat->xstride,
                          bb_ptr);

        if (gmx_debug_at)
        {
            int bbo;
            bbo = (a0 - grid->cell0*grid->na_sc)/grid->na_c;
            fprintf(debug,"%2d %2d %2d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    sx,sy,sz,
                    (grid->bb+bbo*NNBSBB_B)[BBL_X],
                    (grid->bb+bbo*NNBSBB_B)[BBU_X],
                    (grid->bb+bbo*NNBSBB_B)[BBL_Y],
                    (grid->bb+bbo*NNBSBB_B)[BBU_Y],
                    (grid->bb+bbo*NNBSBB_B)[BBL_Z],
                    (grid->bb+bbo*NNBSBB_B)[BBU_Z]);
        }
    }
}

/* Spatially sort the atoms within one grid column */
static void sort_columns_simple(const nbnxn_search_t nbs,
                                int dd_zone,
                                nbnxn_grid_t *grid,
                                int a0,int a1,
                                const int *atinfo,
                                rvec *x,
                                nbnxn_atomdata_t *nbat,
                                int cxy_start,int cxy_end,
                                int *sort_work)
{
    int  cxy;
    int  cx,cy,cz,ncz,cfilled,c;
    int  na,ash,ind,a;
    int  na_c,ash_c;

    if (debug)
    {
        fprintf(debug,"cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0,cxy_start,cxy_end,a0,a1);
    }

    /* Sort the atoms within each x,y column in 3 dimensions */
    for(cxy=cxy_start; cxy<cxy_end; cxy++)
    {
        cx = cxy/grid->ncy;
        cy = cxy - cx*grid->ncy;

        na  = grid->cxy_na[cxy];
        ncz = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,FALSE,
                   nbs->a+ash,na,x,
                   grid->c0[ZZ],
                   ncz*grid->na_sc*SORT_GRID_OVERSIZE/nbs->box[ZZ][ZZ],
                   ncz*grid->na_sc*SGSF,sort_work);

        /* Fill the ncz cells in this column */
        cfilled = grid->cxy_ind[cxy];
        for(cz=0; cz<ncz; cz++)
        {
            c  = grid->cxy_ind[cxy] + cz ;

            ash_c = ash + cz*grid->na_sc;
            na_c  = min(grid->na_sc,na-(ash_c-ash));

            fill_cell(nbs,grid,nbat,
                      ash_c,ash_c+na_c,atinfo,x,
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
            grid->bbcz[c*NNBSBB_D  ] = grid->bb[cfilled*NNBSBB_B+2];
            grid->bbcz[c*NNBSBB_D+1] = grid->bb[cfilled*NNBSBB_B+6];
        }

        /* Set the unused atom indices to -1 */
        for(ind=na; ind<ncz*grid->na_sc; ind++)
        {
            nbs->a[ash+ind] = -1;
        }
    }
}

/* Spatially sort the atoms within one grid column */
static void sort_columns_supersub(const nbnxn_search_t nbs,
                                  int dd_zone,
                                  nbnxn_grid_t *grid,
                                  int a0,int a1,
                                  const int *atinfo,
                                  rvec *x,
                                  nbnxn_atomdata_t *nbat,
                                  int cxy_start,int cxy_end,
                                  int *sort_work)
{
    int  cxy;
    int  cx,cy,cz=-1,c=-1,ncz;
    int  na,ash,na_c,ind,a;
    int  subdiv_z,sub_z,na_z,ash_z;
    int  subdiv_y,sub_y,na_y,ash_y;
    int  subdiv_x,sub_x,na_x,ash_x;

    /* cppcheck-suppress unassignedVariable */
    float bb_work_array[NNBSBB_B+3],*bb_work_align;

    bb_work_align = (float *)(((size_t)(bb_work_array+3)) & (~((size_t)15)));

    if (debug)
    {
        fprintf(debug,"cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0,cxy_start,cxy_end,a0,a1);
    }

    subdiv_x = grid->na_c;
    subdiv_y = GPU_NSUBCELL_X*subdiv_x;
    subdiv_z = GPU_NSUBCELL_Y*subdiv_y;

    /* Sort the atoms within each x,y column in 3 dimensions */
    for(cxy=cxy_start; cxy<cxy_end; cxy++)
    {
        cx = cxy/grid->ncy;
        cy = cxy - cx*grid->ncy;

        na  = grid->cxy_na[cxy];
        ncz = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,FALSE,
                   nbs->a+ash,na,x,
                   grid->c0[ZZ],
                   ncz*grid->na_sc*SORT_GRID_OVERSIZE/nbs->box[ZZ][ZZ],
                   ncz*grid->na_sc*SGSF,sort_work);

        /* This loop goes over the supercells and subcells along z at once */
        for(sub_z=0; sub_z<ncz*GPU_NSUBCELL_Z; sub_z++)
        {
            ash_z = ash + sub_z*subdiv_z;
            na_z  = min(subdiv_z,na-(ash_z-ash));

            /* We have already sorted on z */

            if (sub_z % GPU_NSUBCELL_Z == 0)
            {
                cz = sub_z/GPU_NSUBCELL_Z;
                c  = grid->cxy_ind[cxy] + cz ;

                /* The number of atoms in this supercell */
                na_c = min(grid->na_sc,na-(ash_z-ash));

                grid->nsubc[c] = min(GPU_NSUBCELL,(na_c+grid->na_c-1)/grid->na_c);

                /* Store the z-boundaries of the super cell */
                grid->bbcz[c*NNBSBB_D  ] = x[nbs->a[ash_z]][ZZ];
                grid->bbcz[c*NNBSBB_D+1] = x[nbs->a[ash_z+na_c-1]][ZZ];
            }

#if GPU_NSUBCELL_Y > 1
            /* Sort the atoms along y */
            sort_atoms(YY,(sub_z & 1),
                       nbs->a+ash_z,na_z,x,
                       grid->c0[YY]+cy*grid->sy,grid->inv_sy,
                       subdiv_y*SGSF,sort_work);
#endif

            for(sub_y=0; sub_y<GPU_NSUBCELL_Y; sub_y++)
            {
                ash_y = ash_z + sub_y*subdiv_y;
                na_y  = min(subdiv_y,na-(ash_y-ash));

#if GPU_NSUBCELL_X > 1
                /* Sort the atoms along x */
                sort_atoms(XX,((cz*GPU_NSUBCELL_Y + sub_y) & 1),
                           nbs->a+ash_y,na_y,x,
                           grid->c0[XX]+cx*grid->sx,grid->inv_sx,
                           subdiv_x*SGSF,sort_work);
#endif

                for(sub_x=0; sub_x<GPU_NSUBCELL_X; sub_x++)
                {
                    ash_x = ash_y + sub_x*subdiv_x;
                    na_x  = min(subdiv_x,na-(ash_x-ash));

                    fill_cell(nbs,grid,nbat,
                              ash_x,ash_x+na_x,atinfo,x,
                              grid->na_c*(cx*GPU_NSUBCELL_X+sub_x) + (dd_zone >> 2),
                              grid->na_c*(cy*GPU_NSUBCELL_Y+sub_y) + (dd_zone & 3),
                              grid->na_c*sub_z,
                              bb_work_align);
                }
            }
        }

        /* Set the unused atom indices to -1 */
        for(ind=na; ind<ncz*grid->na_sc; ind++)
        {
            nbs->a[ash+ind] = -1;
        }
    }
}

/* Determine in which grid column atoms should go */
static void calc_column_indices(nbnxn_grid_t *grid,
                                int a0,int a1,
                                rvec *x,const int *move,
                                int thread,int nthread,
                                int *cell,
                                int *cxy_na)
{
    int  n0,n1,i;
    int  cx,cy;

    /* We add one extra cell for particles which moved during DD */
    for(i=0; i<grid->ncx*grid->ncy+1; i++)
    {
        cxy_na[i] = 0;
    }

    n0 = a0 + (int)((thread+0)*(a1 - a0))/nthread;
    n1 = a0 + (int)((thread+1)*(a1 - a0))/nthread;
    for(i=n0; i<n1; i++)
    {
        if (move == NULL || move[i] >= 0)
        {
            /* We need to be careful with rounding,
             * particles might be a few bits outside the local box.
             * The int cast takes care of the lower bound,
             * we need to explicitly take care of the upper bound.
             */
            cx = (int)((x[i][XX] - grid->c0[XX])*grid->inv_sx);
            if (cx == grid->ncx)
            {
                cx = grid->ncx - 1;
            }
            cy = (int)((x[i][YY] - grid->c0[YY])*grid->inv_sy);
            if (cy == grid->ncy)
            {
                cy = grid->ncy - 1;
            }
            /* For the moment cell contains only the, grid local,
             * x and y indices, not z.
             */
            cell[i] = cx*grid->ncy + cy;

#ifdef DEBUG_NBNXN_GRIDDING
            if (cell[i] < 0 || cell[i] >= grid->ncx*grid->ncy)
            {
                gmx_fatal(FARGS,
                          "grid cell cx %d cy %d out of range (max %d %d)",
                          cx,cy,grid->ncx,grid->ncy);
            }
#endif
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

/* Determine in which grid cells the atoms should go */
static void calc_cell_indices(const nbnxn_search_t nbs,
                              int dd_zone,
                              nbnxn_grid_t *grid,
                              int a0,int a1,
                              const int *atinfo,
                              rvec *x,
                              const int *move,
                              nbnxn_atomdata_t *nbat)
{
    int  n0,n1,i;
    int  cx,cy,cxy,ncz_max,ncz;
    int  nthread,thread;
    int  *cxy_na,cxy_na_i;

    nthread = gmx_omp_nthreads_get(emntPairsearch);

#pragma omp parallel for num_threads(nthread) schedule(static)
    for(thread=0; thread<nthread; thread++)
    {
        calc_column_indices(grid,a0,a1,x,move,thread,nthread,
                            nbs->cell,nbs->work[thread].cxy_na);
    }

    /* Make the cell index as a function of x and y */
    ncz_max = 0;
    ncz = 0;
    grid->cxy_ind[0] = 0;
    for(i=0; i<grid->ncx*grid->ncy+1; i++)
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
        for(thread=1; thread<nthread; thread++)
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
        fprintf(debug,"ns na_sc %d na_c %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                grid->na_sc,grid->na_c,grid->nc,
                grid->ncx,grid->ncy,grid->nc/((double)(grid->ncx*grid->ncy)),
                ncz_max);
        if (gmx_debug_at)
        {
            i = 0;
            for(cy=0; cy<grid->ncy; cy++)
            {
                for(cx=0; cx<grid->ncx; cx++)
                {
                    fprintf(debug," %2d",grid->cxy_ind[i+1]-grid->cxy_ind[i]);
                    i++;
                }
                fprintf(debug,"\n");
            }
        }
    }

    /* Make sure the work array for sorting is large enough */
    if (ncz_max*grid->na_sc*SGSF > nbs->work[0].sort_work_nalloc)
    {
        for(thread=0; thread<nbs->nthread_max; thread++)
        {
            nbs->work[thread].sort_work_nalloc =
                over_alloc_large(ncz_max*grid->na_sc*SGSF);
            srenew(nbs->work[thread].sort_work,
                   nbs->work[thread].sort_work_nalloc);
        }
    }

    /* Now we know the dimensions we can fill the grid.
     * This is the first, unsorted fill. We sort the columns after this.
     */
    for(i=a0; i<a1; i++)
    {
        /* At this point nbs->cell contains the local grid x,y indices */
        cxy = nbs->cell[i];
        nbs->a[(grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc + grid->cxy_na[cxy]++] = i;
    }

    /* Set the cell indices for the moved particles */
    n0 = grid->nc*grid->na_sc;
    n1 = grid->nc*grid->na_sc+grid->cxy_na[grid->ncx*grid->ncy];
    for(i=n0; i<n1; i++)
    {
        nbs->cell[nbs->a[i]] = i;
    }

    /* Sort the super-cell columns along z into the sub-cells. */
#pragma omp parallel for num_threads(nbs->nthread_max) schedule(static)
    for(thread=0; thread<nbs->nthread_max; thread++)
    {
        if (grid->bSimple)
        {
            sort_columns_simple(nbs,dd_zone,grid,a0,a1,atinfo,x,nbat,
                                ((thread+0)*grid->ncx*grid->ncy)/nthread,
                                ((thread+1)*grid->ncx*grid->ncy)/nthread,
                                nbs->work[thread].sort_work);
        }
        else
        {
            sort_columns_supersub(nbs,dd_zone,grid,a0,a1,atinfo,x,nbat,
                                  ((thread+0)*grid->ncx*grid->ncy)/nthread,
                                  ((thread+1)*grid->ncx*grid->ncy)/nthread,
                                  nbs->work[thread].sort_work);
        }
    }

#ifdef NBNXN_SEARCH_SSE
    if (grid->bSimple && nbat->XFormat == nbatX8)
    {
        combine_bounding_box_pairs(grid,grid->bb);
    }
#endif

    if (!grid->bSimple)
    {
        grid->nsubc_tot = 0;
        for(i=0; i<grid->nc; i++)
        {
            grid->nsubc_tot += grid->nsubc[i];
        }
    }

    if (debug)
    {
        if (grid->bSimple)
        {
            print_bbsizes_simple(debug,nbs,grid);
        }
        else
        {
            fprintf(debug,"ns non-zero sub-cells: %d average atoms %.2f\n",
                    grid->nsubc_tot,(a1-a0)/(double)grid->nsubc_tot);

            print_bbsizes_supersub(debug,nbs,grid);
        }
    }
}

static void init_grid_flags(nbnxn_cellblock_flags *flags,
                            const nbnxn_grid_t *grid)
{
    int cb;

    flags->ncb = (grid->nc + NBNXN_CELLBLOCK_SIZE - 1)/NBNXN_CELLBLOCK_SIZE;
    if (flags->ncb > flags->flag_nalloc)
    {
        flags->flag_nalloc = over_alloc_large(flags->ncb);
        srenew(flags->flag,flags->flag_nalloc);
    }
    for(cb=0; cb<flags->ncb; cb++)
    {
        flags->flag[cb] = 0;
    }

    flags->bUse = TRUE;
}

/* Sets up a grid and puts the atoms on the grid.
 * This function only operates on one domain of the domain decompostion.
 * Note that without domain decomposition there is only one domain.
 */
void nbnxn_put_on_grid(nbnxn_search_t nbs,
                       int ePBC,matrix box,
                       int dd_zone,
                       rvec corner0,rvec corner1,
                       int a0,int a1,
                       real atom_density,
                       const int *atinfo,
                       rvec *x,
                       int nmoved,int *move,
                       int nb_kernel_type,
                       nbnxn_atomdata_t *nbat)
{
    nbnxn_grid_t *grid;
    int n;
    int nc_max_grid,nc_max;

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
        copy_mat(box,nbs->box);

        if (atom_density >= 0)
        {
            grid->atom_density = atom_density;
        }
        else
        {
            grid->atom_density = grid_atom_density(n-nmoved,corner0,corner1);
        }

        grid->cell0 = 0;

        nbs->natoms_local    = a1 - nmoved;
        /* We assume that nbnxn_put_on_grid is called first
         * for the local atoms (dd_zone=0).
         */
        nbs->natoms_nonlocal = a1 - nmoved;
    }
    else
    {
        nbs->natoms_nonlocal = max(nbs->natoms_nonlocal,a1);
    }

    nc_max_grid = set_grid_size_xy(nbs,grid,n-nmoved,corner0,corner1,
                                   nbs->grid[0].atom_density,
                                   nbat->XFormat);

    nc_max = grid->cell0 + nc_max_grid;

    if (a1 > nbs->cell_nalloc)
    {
        nbs->cell_nalloc = over_alloc_large(a1);
        srenew(nbs->cell,nbs->cell_nalloc);
    }

    /* To avoid conditionals we store the moved particles at the end of a,
     * make sure we have enough space.
     */
    if (nc_max*grid->na_sc + nmoved > nbs->a_nalloc)
    {
        nbs->a_nalloc = over_alloc_large(nc_max*grid->na_sc + nmoved);
        srenew(nbs->a,nbs->a_nalloc);
    }

    if (nc_max*grid->na_sc > nbat->nalloc)
    {
        nbnxn_atomdata_realloc(nbat,nc_max*grid->na_sc);
    }

    calc_cell_indices(nbs,dd_zone,grid,a0,a1,atinfo,x,move,nbat);

    if (dd_zone == 0)
    {
        nbat->natoms_local = nbat->natoms;
    }

    init_grid_flags(&grid->cellblock_flags,grid);

    nbs_cycle_stop(&nbs->cc[enbsCCgrid]);
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nbnxn_search_t nbs,
                                const gmx_domdec_zones_t *zones,
                                const int *atinfo,
                                rvec *x,
                                int nb_kernel_type,
                                nbnxn_atomdata_t *nbat)
{
    int  zone,d;
    rvec c0,c1;

    for(zone=1; zone<zones->n; zone++)
    {
        for(d=0; d<DIM; d++)
        {
            c0[d] = zones->size[zone].bb_x0[d];
            c1[d] = zones->size[zone].bb_x1[d];
        }

        nbnxn_put_on_grid(nbs,nbs->ePBC,NULL,
                          zone,c0,c1,
                          zones->cg_range[zone],
                          zones->cg_range[zone+1],
                          -1,
                          atinfo,
                          x,
                          0,NULL,
                          nb_kernel_type,
                          nbat);
    }
}

/* Add simple grid type information to the local super/sub grid */
void nbnxn_grid_add_simple(nbnxn_search_t nbs,
                           nbnxn_atomdata_t *nbat)
{
    nbnxn_grid_t *grid;
    float *bbcz,*bb;
    int ncd,sc;

    grid = &nbs->grid[0];

    if (grid->bSimple)
    {
        gmx_incons("nbnxn_grid_simple called with a simple grid");
    }

    ncd = grid->na_sc/NBNXN_CPU_CLUSTER_I_SIZE;

    if (grid->nc*ncd > grid->nc_nalloc_simple)
    {
        grid->nc_nalloc_simple = over_alloc_large(grid->nc*ncd);
        srenew(grid->bbcz_simple,grid->nc_nalloc_simple*NNBSBB_D);
        srenew(grid->bb_simple,grid->nc_nalloc_simple*NNBSBB_B);
        srenew(grid->flags_simple,grid->nc_nalloc_simple);
        if (nbat->XFormat)
        {
            sfree_aligned(grid->bbj);
            snew_aligned(grid->bbj,grid->nc_nalloc_simple/2,16);
        }
    }

    bbcz = grid->bbcz_simple;
    bb   = grid->bb_simple;

#pragma omp parallel for num_threads(gmx_omp_nthreads_get(emntPairsearch)) schedule(static)
    for(sc=0; sc<grid->nc; sc++)
    {
        int c,tx,na;

        for(c=0; c<ncd; c++)
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
                    calc_bounding_box_x_x4(na,nbat->x+tx*STRIDE_P4,
                                           bb+tx*NNBSBB_B);
                    break;
                case nbatX8:
                    /* PACK_X8>NBNXN_CPU_CLUSTER_I_SIZE, more complicated */
                    calc_bounding_box_x_x8(na,nbat->x+X8_IND_A(tx*NBNXN_CPU_CLUSTER_I_SIZE),
                                           bb+tx*NNBSBB_B);
                    break;
                default:
                    calc_bounding_box(na,nbat->xstride,
                                      nbat->x+tx*NBNXN_CPU_CLUSTER_I_SIZE*nbat->xstride,
                                      bb+tx*NNBSBB_B);
                    break;
                }
                bbcz[tx*NNBSBB_D+0] = bb[tx*NNBSBB_B         +ZZ];
                bbcz[tx*NNBSBB_D+1] = bb[tx*NNBSBB_B+NNBSBB_C+ZZ];

                /* No interaction optimization yet here */
                grid->flags_simple[tx] = NBNXN_CI_DO_LJ(0) | NBNXN_CI_DO_COUL(0);
            }
            else
            {
                grid->flags_simple[tx] = 0;
            }
        }
    }

#ifdef NBNXN_SEARCH_SSE
    if (grid->bSimple && nbat->XFormat == nbatX8)
    {
        combine_bounding_box_pairs(grid,grid->bb_simple);
    }
#endif
}

void nbnxn_get_ncells(nbnxn_search_t nbs,int *ncx,int *ncy)
{
    *ncx = nbs->grid[0].ncx;
    *ncy = nbs->grid[0].ncy;
}

void nbnxn_get_atomorder(nbnxn_search_t nbs,int **a,int *n)
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
    int ao,cx,cy,cxy,cz,j;

    /* Set the atom order for the home cell (index 0) */
    grid = &nbs->grid[0];

    ao = 0;
    for(cx=0; cx<grid->ncx; cx++)
    {
        for(cy=0; cy<grid->ncy; cy++)
        {
            cxy = cx*grid->ncy + cy;
            j   = grid->cxy_ind[cxy]*grid->na_sc;
            for(cz=0; cz<grid->cxy_na[cxy]; cz++)
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
static void get_cell_range(real b0,real b1,
                           int nc,real c0,real s,real invs,
                           real d2,real r2,int *cf,int *cl)
{
    *cf = max((int)((b0 - c0)*invs),0);

    while (*cf > 0 && d2 + sqr((b0 - c0) - (*cf-1+1)*s) < r2)
    {
        (*cf)--;
    }

    *cl = min((int)((b1 - c0)*invs),nc-1);
    while (*cl < nc-1 && d2 + sqr((*cl+1)*s - (b1 - c0)) < r2)
    {
        (*cl)++;
    }
}

/* Reference code calculating the distance^2 between two bounding boxes */
static float box_dist2(float bx0,float bx1,float by0,
                       float by1,float bz0,float bz1,
                       const float *bb)
{
    float d2;
    float dl,dh,dm,dm0;

    d2 = 0;

    dl  = bx0 - bb[BBU_X];
    dh  = bb[BBL_X] - bx1;
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = by0 - bb[BBU_Y];
    dh  = bb[BBL_Y] - by1;
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = bz0 - bb[BBU_Z];
    dh  = bb[BBL_Z] - bz1;
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    return d2;
}

/* Plain C code calculating the distance^2 between two bounding boxes */
static float subc_bb_dist2(int si,const float *bb_i_ci,
                           int csj,const float *bb_j_all)
{
    const float *bb_i,*bb_j;
    float d2;
    float dl,dh,dm,dm0;

    bb_i = bb_i_ci  +  si*NNBSBB_B;
    bb_j = bb_j_all + csj*NNBSBB_B;

    d2 = 0;

    dl  = bb_i[BBL_X] - bb_j[BBU_X];
    dh  = bb_j[BBL_X] - bb_i[BBU_X];
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = bb_i[BBL_Y] - bb_j[BBU_Y];
    dh  = bb_j[BBL_Y] - bb_i[BBU_Y];
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = bb_i[BBL_Z] - bb_j[BBU_Z];
    dh  = bb_j[BBL_Z] - bb_i[BBU_Z];
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    return d2;
}

#ifdef NBNXN_SEARCH_SSE

/* SSE code for bb distance for bb format xyz0 */
static float subc_bb_dist2_sse(int na_c,
                              int si,const float *bb_i_ci,
                              int csj,const float *bb_j_all)
{
    const float *bb_i,*bb_j;

    __m128 bb_i_SSE0,bb_i_SSE1;
    __m128 bb_j_SSE0,bb_j_SSE1;
    __m128 dl_SSE;
    __m128 dh_SSE;
    __m128 dm_SSE;
    __m128 dm0_SSE;
    __m128 d2_SSE;
#ifndef GMX_X86_SSE4_1
    float d2_array[7],*d2_align;

    d2_align = (float *)(((size_t)(d2_array+3)) & (~((size_t)15)));
#else
    float d2;
#endif

    bb_i = bb_i_ci  +  si*NNBSBB_B;
    bb_j = bb_j_all + csj*NNBSBB_B;

    bb_i_SSE0 = _mm_load_ps(bb_i);
    bb_i_SSE1 = _mm_load_ps(bb_i+NNBSBB_C);
    bb_j_SSE0 = _mm_load_ps(bb_j);
    bb_j_SSE1 = _mm_load_ps(bb_j+NNBSBB_C);

    dl_SSE    = _mm_sub_ps(bb_i_SSE0,bb_j_SSE1);
    dh_SSE    = _mm_sub_ps(bb_j_SSE0,bb_i_SSE1);

    dm_SSE    = _mm_max_ps(dl_SSE,dh_SSE);
    dm0_SSE   = _mm_max_ps(dm_SSE,_mm_setzero_ps());
#ifndef GMX_X86_SSE4_1
    d2_SSE    = _mm_mul_ps(dm0_SSE,dm0_SSE);

    _mm_store_ps(d2_align,d2_SSE);

    return d2_align[0] + d2_align[1] + d2_align[2];
#else
    /* SSE4.1 dot product of components 0,1,2 */
    d2_SSE    = _mm_dp_ps(dm0_SSE,dm0_SSE,0x71);

    _mm_store_ss(&d2,d2_SSE);

    return d2;
#endif
}

/* Calculate bb bounding distances of bb_i[si,...,si+3] and store them in d2 */
#define SUBC_BB_DIST2_SSE_XXXX_INNER(si,bb_i,d2) \
{                                                \
    int    shi;                                  \
                                                 \
    __m128 dx_0,dy_0,dz_0;                       \
    __m128 dx_1,dy_1,dz_1;                       \
                                                 \
    __m128 mx,my,mz;                             \
    __m128 m0x,m0y,m0z;                          \
                                                 \
    __m128 d2x,d2y,d2z;                          \
    __m128 d2s,d2t;                              \
                                                 \
    shi = si*NNBSBB_D*DIM;                       \
                                                 \
    xi_l = _mm_load_ps(bb_i+shi+0*STRIDE_8BB);   \
    yi_l = _mm_load_ps(bb_i+shi+1*STRIDE_8BB);   \
    zi_l = _mm_load_ps(bb_i+shi+2*STRIDE_8BB);   \
    xi_h = _mm_load_ps(bb_i+shi+3*STRIDE_8BB);   \
    yi_h = _mm_load_ps(bb_i+shi+4*STRIDE_8BB);   \
    zi_h = _mm_load_ps(bb_i+shi+5*STRIDE_8BB);   \
                                                 \
    dx_0 = _mm_sub_ps(xi_l,xj_h);                \
    dy_0 = _mm_sub_ps(yi_l,yj_h);                \
    dz_0 = _mm_sub_ps(zi_l,zj_h);                \
                                                 \
    dx_1 = _mm_sub_ps(xj_l,xi_h);                \
    dy_1 = _mm_sub_ps(yj_l,yi_h);                \
    dz_1 = _mm_sub_ps(zj_l,zi_h);                \
                                                 \
    mx   = _mm_max_ps(dx_0,dx_1);                \
    my   = _mm_max_ps(dy_0,dy_1);                \
    mz   = _mm_max_ps(dz_0,dz_1);                \
                                                 \
    m0x  = _mm_max_ps(mx,zero);                  \
    m0y  = _mm_max_ps(my,zero);                  \
    m0z  = _mm_max_ps(mz,zero);                  \
                                                 \
    d2x  = _mm_mul_ps(m0x,m0x);                  \
    d2y  = _mm_mul_ps(m0y,m0y);                  \
    d2z  = _mm_mul_ps(m0z,m0z);                  \
                                                 \
    d2s  = _mm_add_ps(d2x,d2y);                  \
    d2t  = _mm_add_ps(d2s,d2z);                  \
                                                 \
    _mm_store_ps(d2+si,d2t);                     \
}

/* SSE code for nsi bb distances for bb format xxxxyyyyzzzz */
static void subc_bb_dist2_sse_xxxx(const float *bb_j,
                                   int nsi,const float *bb_i,
                                   float *d2)
{
    __m128 xj_l,yj_l,zj_l;
    __m128 xj_h,yj_h,zj_h;
    __m128 xi_l,yi_l,zi_l;
    __m128 xi_h,yi_h,zi_h;

    __m128 zero;

    zero = _mm_setzero_ps();

    xj_l = _mm_set1_ps(bb_j[0*STRIDE_8BB]);
    yj_l = _mm_set1_ps(bb_j[1*STRIDE_8BB]);
    zj_l = _mm_set1_ps(bb_j[2*STRIDE_8BB]);
    xj_h = _mm_set1_ps(bb_j[3*STRIDE_8BB]);
    yj_h = _mm_set1_ps(bb_j[4*STRIDE_8BB]);
    zj_h = _mm_set1_ps(bb_j[5*STRIDE_8BB]);

    /* Here we "loop" over si (0,STRIDE_8BB) from 0 to nsi with step STRIDE_8BB.
     * But as we know the number of iterations is 1 or 2, we unroll manually.
     */
    SUBC_BB_DIST2_SSE_XXXX_INNER(0,bb_i,d2);
    if (STRIDE_8BB < nsi)
    {
        SUBC_BB_DIST2_SSE_XXXX_INNER(STRIDE_8BB,bb_i,d2);
    }
}

#endif /* NBNXN_SEARCH_SSE */

/* Plain C function which determines if any atom pair between two cells
 * is within distance sqrt(rl2).
 */
static gmx_bool subc_in_range_x(int na_c,
                                int si,const real *x_i,
                                int csj,int stride,const real *x_j,
                                real rl2)
{
    int  i,j,i0,j0;
    real d2;

    for(i=0; i<na_c; i++)
    {
        i0 = (si*na_c + i)*DIM;
        for(j=0; j<na_c; j++)
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

/* SSE function which determines if any atom pair between two cells,
 * both with 8 atoms, is within distance sqrt(rl2).
 */
static gmx_bool subc_in_range_sse8(int na_c,
                                   int si,const real *x_i,
                                   int csj,int stride,const real *x_j,
                                   real rl2)
{
#ifdef NBNXN_SEARCH_SSE_SINGLE
    __m128 ix_SSE0,iy_SSE0,iz_SSE0;
    __m128 ix_SSE1,iy_SSE1,iz_SSE1;

    __m128 rc2_SSE;

    int na_c_sse;
    int j0,j1;

    rc2_SSE   = _mm_set1_ps(rl2);

    na_c_sse = NBNXN_GPU_CLUSTER_SIZE/STRIDE_8BB;
    ix_SSE0 = _mm_load_ps(x_i+(si*na_c_sse*DIM+0)*STRIDE_8BB);
    iy_SSE0 = _mm_load_ps(x_i+(si*na_c_sse*DIM+1)*STRIDE_8BB);
    iz_SSE0 = _mm_load_ps(x_i+(si*na_c_sse*DIM+2)*STRIDE_8BB);
    ix_SSE1 = _mm_load_ps(x_i+(si*na_c_sse*DIM+3)*STRIDE_8BB);
    iy_SSE1 = _mm_load_ps(x_i+(si*na_c_sse*DIM+4)*STRIDE_8BB);
    iz_SSE1 = _mm_load_ps(x_i+(si*na_c_sse*DIM+5)*STRIDE_8BB);

    /* We loop from the outer to the inner particles to maximize
     * the chance that we find a pair in range quickly and return.
     */
    j0 = csj*na_c;
    j1 = j0 + na_c - 1;
    while (j0 < j1)
    {
        __m128 jx0_SSE,jy0_SSE,jz0_SSE;
        __m128 jx1_SSE,jy1_SSE,jz1_SSE;

        __m128 dx_SSE0,dy_SSE0,dz_SSE0;
        __m128 dx_SSE1,dy_SSE1,dz_SSE1;
        __m128 dx_SSE2,dy_SSE2,dz_SSE2;
        __m128 dx_SSE3,dy_SSE3,dz_SSE3;

        __m128 rsq_SSE0;
        __m128 rsq_SSE1;
        __m128 rsq_SSE2;
        __m128 rsq_SSE3;

        __m128 wco_SSE0;
        __m128 wco_SSE1;
        __m128 wco_SSE2;
        __m128 wco_SSE3;
        __m128 wco_any_SSE01,wco_any_SSE23,wco_any_SSE;

        jx0_SSE = _mm_load1_ps(x_j+j0*stride+0);
        jy0_SSE = _mm_load1_ps(x_j+j0*stride+1);
        jz0_SSE = _mm_load1_ps(x_j+j0*stride+2);

        jx1_SSE = _mm_load1_ps(x_j+j1*stride+0);
        jy1_SSE = _mm_load1_ps(x_j+j1*stride+1);
        jz1_SSE = _mm_load1_ps(x_j+j1*stride+2);

        /* Calculate distance */
        dx_SSE0            = _mm_sub_ps(ix_SSE0,jx0_SSE);
        dy_SSE0            = _mm_sub_ps(iy_SSE0,jy0_SSE);
        dz_SSE0            = _mm_sub_ps(iz_SSE0,jz0_SSE);
        dx_SSE1            = _mm_sub_ps(ix_SSE1,jx0_SSE);
        dy_SSE1            = _mm_sub_ps(iy_SSE1,jy0_SSE);
        dz_SSE1            = _mm_sub_ps(iz_SSE1,jz0_SSE);
        dx_SSE2            = _mm_sub_ps(ix_SSE0,jx1_SSE);
        dy_SSE2            = _mm_sub_ps(iy_SSE0,jy1_SSE);
        dz_SSE2            = _mm_sub_ps(iz_SSE0,jz1_SSE);
        dx_SSE3            = _mm_sub_ps(ix_SSE1,jx1_SSE);
        dy_SSE3            = _mm_sub_ps(iy_SSE1,jy1_SSE);
        dz_SSE3            = _mm_sub_ps(iz_SSE1,jz1_SSE);

        /* rsq = dx*dx+dy*dy+dz*dz */
        rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0,dy_SSE0,dz_SSE0);
        rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1,dy_SSE1,dz_SSE1);
        rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2,dy_SSE2,dz_SSE2);
        rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3,dy_SSE3,dz_SSE3);

        wco_SSE0           = _mm_cmplt_ps(rsq_SSE0,rc2_SSE);
        wco_SSE1           = _mm_cmplt_ps(rsq_SSE1,rc2_SSE);
        wco_SSE2           = _mm_cmplt_ps(rsq_SSE2,rc2_SSE);
        wco_SSE3           = _mm_cmplt_ps(rsq_SSE3,rc2_SSE);

        wco_any_SSE01      = _mm_or_ps(wco_SSE0,wco_SSE1);
        wco_any_SSE23      = _mm_or_ps(wco_SSE2,wco_SSE3);
        wco_any_SSE        = _mm_or_ps(wco_any_SSE01,wco_any_SSE23);

        if (_mm_movemask_ps(wco_any_SSE))
        {
            return TRUE;
        }

        j0++;
        j1--;
    }
    return FALSE;

#else
    /* No SSE */
    gmx_incons("SSE function called without SSE support");

    return TRUE;
#endif
}

/* Returns the j sub-cell for index cj_ind */
static int nbl_cj(const nbnxn_pairlist_t *nbl,int cj_ind)
{
    return nbl->cj4[cj_ind>>2].cj[cj_ind & 3];
}

/* Returns the i-interaction mask of the j sub-cell for index cj_ind */
static unsigned nbl_imask0(const nbnxn_pairlist_t *nbl,int cj_ind)
{
    return nbl->cj4[cj_ind>>2].imei[0].imask;
}

/* Ensures there is enough space for extra extra exclusion masks */
static void check_excl_space(nbnxn_pairlist_t *nbl,int extra)
{
    if (nbl->nexcl+extra > nbl->excl_nalloc)
    {
        nbl->excl_nalloc = over_alloc_small(nbl->nexcl+extra);
        nbnxn_realloc_void((void **)&nbl->excl,
                           nbl->nexcl*sizeof(*nbl->excl),
                           nbl->excl_nalloc*sizeof(*nbl->excl),
                           nbl->alloc,nbl->free);
    }
}

/* Ensures there is enough space for ncell extra j-cells in the list */
static void check_subcell_list_space_simple(nbnxn_pairlist_t *nbl,
                                            int ncell)
{
    int cj_max;

    cj_max = nbl->ncj + ncell;

    if (cj_max > nbl->cj_nalloc)
    {
        nbl->cj_nalloc = over_alloc_small(cj_max);
        nbnxn_realloc_void((void **)&nbl->cj,
                           nbl->ncj*sizeof(*nbl->cj),
                           nbl->cj_nalloc*sizeof(*nbl->cj),
                           nbl->alloc,nbl->free);
    }
}

/* Ensures there is enough space for ncell extra j-subcells in the list */
static void check_subcell_list_space_supersub(nbnxn_pairlist_t *nbl,
                                              int nsupercell)
{
    int ncj4_max,j4,j,w,t;

#define NWARP       2
#define WARP_SIZE  32

    /* We can have maximally nsupercell*GPU_NSUBCELL sj lists */
    /* We can store 4 j-subcell - i-supercell pairs in one struct.
     * since we round down, we need one extra entry.
     */
    ncj4_max = ((nbl->work->cj_ind + nsupercell*GPU_NSUBCELL + 4-1) >> 2);

    if (ncj4_max > nbl->cj4_nalloc)
    {
        nbl->cj4_nalloc = over_alloc_small(ncj4_max);
        nbnxn_realloc_void((void **)&nbl->cj4,
                           nbl->work->cj4_init*sizeof(*nbl->cj4),
                           nbl->cj4_nalloc*sizeof(*nbl->cj4),
                           nbl->alloc,nbl->free);
    }

    if (ncj4_max > nbl->work->cj4_init)
    {
        for(j4=nbl->work->cj4_init; j4<ncj4_max; j4++)
        {
            /* No i-subcells and no excl's in the list initially */
            for(w=0; w<NWARP; w++)
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

    for(t=0; t<WARP_SIZE; t++)
    {
        /* Turn all interaction bits on */
        excl->pair[t] = NBNXN_INT_MASK_ALL;
    }
}

/* Initializes a single nbnxn_pairlist_t data structure */
static void nbnxn_init_pairlist(nbnxn_pairlist_t *nbl,
                                gmx_bool bSimple,
                                nbnxn_alloc_t *alloc,
                                nbnxn_free_t  *free)
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
        check_excl_space(nbl,1);
        nbl->nexcl       = 1;
        set_no_excls(&nbl->excl[0]);
    }

    snew(nbl->work,1);
#ifdef NBNXN_BBXXXX
    snew_aligned(nbl->work->bb_ci,GPU_NSUBCELL/STRIDE_8BB*NNBSBB_XXXX,16);
#else
    snew_aligned(nbl->work->bb_ci,GPU_NSUBCELL*NNBSBB_B,16);
#endif
    snew_aligned(nbl->work->x_ci,NBNXN_NA_SC_MAX*DIM,16);
#ifdef NBNXN_SEARCH_SSE
    snew_aligned(nbl->work->x_ci_x86_simd128,1,16);
#ifdef GMX_X86_AVX_256
    snew_aligned(nbl->work->x_ci_x86_simd256,1,32);
#endif
#endif
    snew_aligned(nbl->work->d2,GPU_NSUBCELL,16);
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
        nbl_list->nnbl > NBNXN_CELLBLOCK_MAX_THREADS)
    {
        gmx_fatal(FARGS,"%d OpenMP threads were requested. Since the non-bonded force buffer reduction is prohibitively slow with more than %d threads, we do not allow this. Use %d or less OpenMP threads.",
                  nbl_list->nnbl,NBNXN_CELLBLOCK_MAX_THREADS,NBNXN_CELLBLOCK_MAX_THREADS);
    }

    snew(nbl_list->nbl,nbl_list->nnbl);
    /* Execute in order to avoid memory interleaving between threads */
#pragma omp parallel for num_threads(nbl_list->nnbl) schedule(static)
    for(i=0; i<nbl_list->nnbl; i++)
    {
        /* Allocate the nblist data structure locally on each thread
         * to optimize memory access for NUMA architectures.
         */
        snew(nbl_list->nbl[i],1);

        /* Only list 0 is used on the GPU, use normal allocation for i>0 */
        if (i == 0)
        {
            nbnxn_init_pairlist(nbl_list->nbl[i],nbl_list->bSimple,alloc,free);
        }
        else
        {
            nbnxn_init_pairlist(nbl_list->nbl[i],nbl_list->bSimple,NULL,NULL);
        }
    }
}

/* Print statistics of a pair list, used for debug output */
static void print_nblist_statistics_simple(FILE *fp,const nbnxn_pairlist_t *nbl,
                                           const nbnxn_search_t nbs,real rl)
{
    const nbnxn_grid_t *grid;
    int cs[SHIFTS];
    int s,i,j;
    int npexcl;

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp,"nbl nci %d ncj %d\n",
            nbl->nci,nbl->ncj);
    fprintf(fp,"nbl na_sc %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->na_sc,rl,nbl->ncj,nbl->ncj/(double)grid->nc,
            nbl->ncj/(double)grid->nc*grid->na_sc,
            nbl->ncj/(double)grid->nc*grid->na_sc/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nc*grid->na_sc/det(nbs->box)));

    fprintf(fp,"nbl average j cell list length %.1f\n",
            0.25*nbl->ncj/(double)nbl->nci);

    for(s=0; s<SHIFTS; s++)
    {
        cs[s] = 0;
    }
    npexcl = 0;
    for(i=0; i<nbl->nci; i++)
    {
        cs[nbl->ci[i].shift & NBNXN_CI_SHIFT] +=
            nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start;

        j = nbl->ci[i].cj_ind_start;
        while (j < nbl->ci[i].cj_ind_end &&
               nbl->cj[j].excl != NBNXN_INT_MASK_ALL)
        {
            npexcl++;
            j++;
        }
    }
    fprintf(fp,"nbl cell pairs, total: %d excl: %d %.1f%%\n",
            nbl->ncj,npexcl,100*npexcl/(double)nbl->ncj);
    for(s=0; s<SHIFTS; s++)
    {
        if (cs[s] > 0)
        {
            fprintf(fp,"nbl shift %2d ncj %3d\n",s,cs[s]);
        }
    }
}

/* Print statistics of a pair lists, used for debug output */
static void print_nblist_statistics_supersub(FILE *fp,const nbnxn_pairlist_t *nbl,
                                             const nbnxn_search_t nbs,real rl)
{
    const nbnxn_grid_t *grid;
    int i,j4,j,si,b;
    int c[GPU_NSUBCELL+1];

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp,"nbl nsci %d ncj4 %d nsi %d excl4 %d\n",
            nbl->nsci,nbl->ncj4,nbl->nci_tot,nbl->nexcl);
    fprintf(fp,"nbl na_c %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->na_ci,rl,nbl->nci_tot,nbl->nci_tot/(double)grid->nsubc_tot,
            nbl->nci_tot/(double)grid->nsubc_tot*grid->na_c,
            nbl->nci_tot/(double)grid->nsubc_tot*grid->na_c/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nsubc_tot*grid->na_c/det(nbs->box)));

    fprintf(fp,"nbl average j super cell list length %.1f\n",
            0.25*nbl->ncj4/(double)nbl->nsci);
    fprintf(fp,"nbl average i sub cell list length %.1f\n",
            nbl->nci_tot/(0.25*nbl->ncj4));

    for(si=0; si<=GPU_NSUBCELL; si++)
    {
        c[si] = 0;
    }
    for(i=0; i<nbl->nsci; i++)
    {
        for(j4=nbl->sci[i].cj4_ind_start; j4<nbl->sci[i].cj4_ind_end; j4++)
        {
            for(j=0; j<4; j++)
            {
                b = 0;
                for(si=0; si<GPU_NSUBCELL; si++)
                {
                    if (nbl->cj4[j4].imei[0].imask & (1U << (j*GPU_NSUBCELL + si)))
                    {
                        b++;
                    }
                }
                c[b]++;
            }
        }
    }
    for(b=0; b<=GPU_NSUBCELL; b++)
    {
        fprintf(fp,"nbl j-list #i-subcell %d %7d %4.1f\n",
                b,c[b],100.0*c[b]/(double)(nbl->ncj4*NBNXN_GPU_JGROUP_SIZE));
    }
}

/* Print the full pair list, used for debug output */
static void print_supersub_nsp(const char *fn,
                               const nbnxn_pairlist_t *nbl,
                               int iloc)
{
    char buf[STRLEN];
    FILE *fp;
    int i,nsp,j4,p;

    sprintf(buf,"%s_%s.xvg",fn,NONLOCAL_I(iloc) ? "nl" : "l");
    fp = ffopen(buf,"w");

    for(i=0; i<nbl->nci; i++)
    {
        nsp = 0;
        for(j4=nbl->sci[i].cj4_ind_start; j4<nbl->sci[i].cj4_ind_end; j4++)
        {
            for(p=0; p<NBNXN_GPU_JGROUP_SIZE*GPU_NSUBCELL; p++)
            {
                nsp += (nbl->cj4[j4].imei[0].imask >> p) & 1;
            }
        }
        fprintf(fp,"%4d %3d %3d\n",
                i,
                nsp,
                nbl->sci[i].cj4_ind_end-nbl->sci[i].cj4_ind_start);
    }

    fclose(fp);
}

/* Returns a pointer to the exclusion mask for cj4-unit cj4, warp warp */
static void low_get_nbl_exclusions(nbnxn_pairlist_t *nbl,int cj4,
                                   int warp,nbnxn_excl_t **excl)
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
 * allocates extra memory, if necessary.
 */
static void get_nbl_exclusions_1(nbnxn_pairlist_t *nbl,int cj4,
                                 int warp,nbnxn_excl_t **excl)
{
    if (nbl->cj4[cj4].imei[warp].excl_ind == 0)
    {
        /* We need to make a new list entry, check if we have space */
        check_excl_space(nbl,1);
    }
    low_get_nbl_exclusions(nbl,cj4,warp,excl);
}

/* Returns pointers to the exclusion mask for cj4-unit cj4 for both warps,
 * allocates extra memory, if necessary.
 */
static void get_nbl_exclusions_2(nbnxn_pairlist_t *nbl,int cj4,
                                 nbnxn_excl_t **excl_w0,
                                 nbnxn_excl_t **excl_w1)
{
    /* Check for space we might need */
    check_excl_space(nbl,2);

    low_get_nbl_exclusions(nbl,cj4,0,excl_w0);
    low_get_nbl_exclusions(nbl,cj4,1,excl_w1);
}

/* Sets the self exclusions i=j and pair exclusions i>j */
static void set_self_and_newton_excls_supersub(nbnxn_pairlist_t *nbl,
                                               int cj4_ind,int sj_offset,
                                               int si)
{
    nbnxn_excl_t *excl[2];
    int  ei,ej,w;

    /* Here we only set the set self and double pair exclusions */

    get_nbl_exclusions_2(nbl,cj4_ind,&excl[0],&excl[1]);

    /* Only minor < major bits set */
    for(ej=0; ej<nbl->na_ci; ej++)
    {
        w = (ej>>2);
        for(ei=ej; ei<nbl->na_ci; ei++)
        {
            excl[w]->pair[(ej&(4-1))*nbl->na_ci+ei] &=
                ~(1U << (sj_offset*GPU_NSUBCELL+si));
        }
    }
}

/* Returns a diagonal or off-diagonal interaction mask for plain C lists */
static unsigned int get_imask(gmx_bool rdiag,int ci,int cj)
{
    return (rdiag && ci == cj ? NBNXN_INT_MASK_DIAG : NBNXN_INT_MASK_ALL);
}

#ifdef NBNXN_SEARCH_SSE
/* Returns a diagonal or off-diagonal interaction mask for SIMD128 lists */
static unsigned int get_imask_x86_simd128(gmx_bool rdiag,int ci,int cj)
{
#ifndef GMX_DOUBLE /* cj-size = 4 */
    return (rdiag && ci == cj ? NBNXN_INT_MASK_DIAG : NBNXN_INT_MASK_ALL);
#else              /* cj-size = 2 */
    return (rdiag && ci*2 == cj ? NBNXN_INT_MASK_DIAG_J2_0 :
            (rdiag && ci*2+1 == cj ? NBNXN_INT_MASK_DIAG_J2_1 :
             NBNXN_INT_MASK_ALL));
#endif
}

#ifdef GMX_X86_AVX_256
/* Returns a diagonal or off-diagonal interaction mask for SIMD256 lists */
static unsigned int get_imask_x86_simd256(gmx_bool rdiag,int ci,int cj)
{
#ifndef GMX_DOUBLE /* cj-size = 8 */
    return (rdiag && ci == cj*2 ? NBNXN_INT_MASK_DIAG_J8_0 :
            (rdiag && ci == cj*2+1 ? NBNXN_INT_MASK_DIAG_J8_1 :
             NBNXN_INT_MASK_ALL));
#else              /* cj-size = 2 */
    return (rdiag && ci == cj ? NBNXN_INT_MASK_DIAG : NBNXN_INT_MASK_ALL);
#endif
}
#endif
#endif /* NBNXN_SEARCH_SSE */

/* Plain C code for making a pair list of cell ci vs cell cjf-cjl.
 * Checks bounding box distances and possibly atom pair distances.
 */
static void make_cluster_list_simple(const nbnxn_grid_t *gridj,
                                     nbnxn_pairlist_t *nbl,
                                     int ci,int cjf,int cjl,
                                     gmx_bool remove_sub_diag,
                                     const real *x_j,
                                     real rl2,float rbb2,
                                     int *ndistc)
{
    const nbnxn_list_work_t *work;

    const float *bb_ci;
    const real  *x_ci;

    gmx_bool   InRange;
    real       d2;
    int        cjf_gl,cjl_gl,cj;

    work = nbl->work;

    bb_ci = nbl->work->bb_ci;
    x_ci  = nbl->work->x_ci;

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
        d2 = subc_bb_dist2(0,bb_ci,cjf,gridj->bb);
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
            int i,j;

            cjf_gl = gridj->cell0 + cjf;
            for(i=0; i<NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for(j=0; j<NBNXN_CPU_CLUSTER_I_SIZE; j++)
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
        d2 = subc_bb_dist2(0,bb_ci,cjl,gridj->bb);
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
            int i,j;

            cjl_gl = gridj->cell0 + cjl;
            for(i=0; i<NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for(j=0; j<NBNXN_CPU_CLUSTER_I_SIZE; j++)
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
        for(cj=cjf; cj<=cjl; cj++)
        {
            /* Store cj and the interaction mask */
            nbl->cj[nbl->ncj].cj   = gridj->cell0 + cj;
            nbl->cj[nbl->ncj].excl = get_imask(remove_sub_diag,ci,cj);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

#ifdef NBNXN_SEARCH_SSE
/* Include make_cluster_list_x86_simd128/256 */
#define GMX_MM128_HERE
#include "gmx_x86_simd_macros.h"
#define STRIDE_S  PACK_X4
#include "nbnxn_search_x86_simd.h"
#undef STRIDE_S
#undef GMX_MM128_HERE
#ifdef GMX_X86_AVX_256
/* Include make_cluster_list_x86_simd128/256 */
#define GMX_MM256_HERE
#include "gmx_x86_simd_macros.h"
#define STRIDE_S  GMX_X86_SIMD_WIDTH_HERE
#include "nbnxn_search_x86_simd.h"
#undef STRIDE_S
#undef GMX_MM256_HERE
#endif
#endif

/* Plain C or SSE code for making a pair list of super-cell sci vs scj.
 * Checks bounding box distances and possibly atom pair distances.
 */
static void make_cluster_list_supersub(const nbnxn_search_t nbs,
                                       const nbnxn_grid_t *gridi,
                                       const nbnxn_grid_t *gridj,
                                       nbnxn_pairlist_t *nbl,
                                       int sci,int scj,
                                       gmx_bool sci_equals_scj,
                                       int stride,const real *x,
                                       real rl2,float rbb2,
                                       int *ndistc)
{
    int  na_c;
    int  npair;
    int  cjo,ci1,ci,cj,cj_gl;
    int  cj4_ind,cj_offset;
    unsigned imask;
    nbnxn_cj4_t *cj4;
    const float *bb_ci;
    const real *x_ci;
    float *d2l,d2;
    int  w;
#define PRUNE_LIST_CPU_ONE
#ifdef PRUNE_LIST_CPU_ONE
    int  ci_last=-1;
#endif

    d2l = nbl->work->d2;

    bb_ci = nbl->work->bb_ci;
    x_ci  = nbl->work->x_ci;

    na_c = gridj->na_c;

    for(cjo=0; cjo<gridj->nsubc[scj]; cjo++)
    {
        cj4_ind   = (nbl->work->cj_ind >> 2);
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
        /* Determine all ci1 bb distances in one call with SSE */
        subc_bb_dist2_sse_xxxx(gridj->bb+(cj>>STRIDE_8BB_2LOG)*NNBSBB_XXXX+(cj & (STRIDE_8BB-1)),
                               ci1,bb_ci,d2l);
        *ndistc += na_c*2;
#endif

        npair = 0;
        /* We use a fixed upper-bound instead of ci1 to help optimization */
        for(ci=0; ci<GPU_NSUBCELL; ci++)
        {
            if (ci == ci1)
            {
                break;
            }

#ifndef NBNXN_BBXXXX
            /* Determine the bb distance between ci and cj */
            d2l[ci] = subc_bb_dist2(ci,bb_ci,cj,gridj->bb);
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
                (d2 < rl2 && subc_in_range_x(na_c,ci,x_ci,cj_gl,stride,x,rl2)))
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
#ifdef NBNXN_8BB_SSE
                !subc_in_range_sse8
#else
                !subc_in_range_x
#endif
                                (na_c,ci_last,x_ci,cj_gl,stride,x,rl2))
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
                set_self_and_newton_excls_supersub(nbl,cj4_ind,cj_offset,cjo);
            }

            /* Copy the cluster interaction mask to the list */
            for(w=0; w<NWARP; w++)
            {
                cj4->imei[w].imask |= imask;
            }

            nbl->work->cj_ind++;

            /* Keep the count */
            nbl->nci_tot += npair;

            /* Increase the closing index in i super-cell list */
            nbl->sci[nbl->nsci].cj4_ind_end = ((nbl->work->cj_ind+4-1)>>2);
        }
    }
}

/* Set all atom-pair exclusions from the topology stored in excl
 * as masks in the pair-list for simple list i-entry nbl_ci
 */
static void set_ci_top_excls(const nbnxn_search_t nbs,
                             nbnxn_pairlist_t *nbl,
                             gmx_bool diagRemoved,
                             int na_ci_2log,
                             int na_cj_2log,
                             const nbnxn_ci_t *nbl_ci,
                             const t_blocka *excl)
{
    const int *cell;
    int ci;
    int cj_ind_first,cj_ind_last;
    int cj_first,cj_last;
    int ndirect;
    int i,ai,aj,si,eind,ge,se;
    int found,cj_ind_0,cj_ind_1,cj_ind_m;
    int cj_m;
    gmx_bool Found_si;
    int si_ind;
    nbnxn_excl_t *nbl_excl;
    int inner_i,inner_e;

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
#ifdef NBNXN_SEARCH_SSE
    else
    {
        while (cj_ind_first + ndirect <= cj_ind_last &&
               nbl->cj[cj_ind_first+ndirect].cj == ci_to_cj(na_cj_2log,ci) + ndirect)
        {
            ndirect++;
        }
    }
#endif

    /* Loop over the atoms in the i super-cell */
    for(i=0; i<nbl->na_sc; i++)
    {
        ai = nbs->a[ci*nbl->na_sc+i];
        if (ai >= 0)
        {
            si  = (i>>na_ci_2log);

            /* Loop over the topology-based exclusions for this i-atom */
            for(eind=excl->index[ai]; eind<excl->index[ai+1]; eind++)
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
                        found = -1;
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

/* Set all atom-pair exclusions from the topology stored in excl
 * as masks in the pair-list for i-super-cell entry nbl_sci
 */
static void set_sci_top_excls(const nbnxn_search_t nbs,
                              nbnxn_pairlist_t *nbl,
                              gmx_bool diagRemoved,
                              int na_c_2log,
                              const nbnxn_sci_t *nbl_sci,
                              const t_blocka *excl)
{
    const int *cell;
    int na_c;
    int sci;
    int cj_ind_first,cj_ind_last;
    int cj_first,cj_last;
    int ndirect;
    int i,ai,aj,si,eind,ge,se;
    int found,cj_ind_0,cj_ind_1,cj_ind_m;
    int cj_m;
    gmx_bool Found_si;
    int si_ind;
    nbnxn_excl_t *nbl_excl;
    int inner_i,inner_e,w;

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
    cj_last  = nbl_cj(nbl,cj_ind_last);

    /* Determine how many contiguous j-clusters we have starting
     * from the first i-cluster. This number can be used to directly
     * calculate j-cluster indices for excluded atoms.
     */
    ndirect = 0;
    while (cj_ind_first + ndirect <= cj_ind_last &&
           nbl_cj(nbl,cj_ind_first+ndirect) == sci*GPU_NSUBCELL + ndirect)
    {
        ndirect++;
    }

    /* Loop over the atoms in the i super-cell */
    for(i=0; i<nbl->na_sc; i++)
    {
        ai = nbs->a[sci*nbl->na_sc+i];
        if (ai >= 0)
        {
            si  = (i>>na_c_2log);

            /* Loop over the topology-based exclusions for this i-atom */
            for(eind=excl->index[ai]; eind<excl->index[ai+1]; eind++)
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
                        found = -1;
                        cj_ind_0 = cj_ind_first + ndirect;
                        cj_ind_1 = cj_ind_last + 1;
                        while (found == -1 && cj_ind_0 < cj_ind_1)
                        {
                            cj_ind_m = (cj_ind_0 + cj_ind_1)>>1;

                            cj_m = nbl_cj(nbl,cj_ind_m);

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

/* Macro for getting the index of atom a within a cluster */
#define AMODI(a)  ((a) & (NBNXN_CPU_CLUSTER_I_SIZE - 1))
/* Macro for converting an atom number to a cluster number */
#define A2CI(a)   ((a) >> NBNXN_CPU_CLUSTER_I_SIZE_2LOG)

                        if (nbl_imask0(nbl,found) & (1U << (AMODI(found)*GPU_NSUBCELL + si)))
                        {
                            w       = (inner_e >> 2);

                            get_nbl_exclusions_1(nbl,A2CI(found),w,&nbl_excl);

                            nbl_excl->pair[AMODI(inner_e)*nbl->na_ci+inner_i] &=
                                ~(1U << (AMODI(found)*GPU_NSUBCELL + si));
                        }

#undef AMODI
#undef A2CI
                    }
                }
            }
        }
    }
}

/* Reallocate the simple ci list for at least n entries */
static void nb_realloc_ci(nbnxn_pairlist_t *nbl,int n)
{
    nbl->ci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->ci,
                       nbl->nci*sizeof(*nbl->ci),
                       nbl->ci_nalloc*sizeof(*nbl->ci),
                       nbl->alloc,nbl->free);
}

/* Reallocate the super-cell sci list for at least n entries */
static void nb_realloc_sci(nbnxn_pairlist_t *nbl,int n)
{
    nbl->sci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->sci,
                       nbl->nsci*sizeof(*nbl->sci),
                       nbl->sci_nalloc*sizeof(*nbl->sci),
                       nbl->alloc,nbl->free);
}

/* Make a new ci entry at index nbl->nci */
static void new_ci_entry(nbnxn_pairlist_t *nbl,int ci,int shift,int flags,
                         nbnxn_list_work_t *work)
{
    if (nbl->nci + 1 > nbl->ci_nalloc)
    {
        nb_realloc_ci(nbl,nbl->nci+1);
    }
    nbl->ci[nbl->nci].ci            = ci;
    nbl->ci[nbl->nci].shift         = shift;
    /* Store the interaction flags along with the shift */
    nbl->ci[nbl->nci].shift        |= flags;
    nbl->ci[nbl->nci].cj_ind_start  = nbl->ncj;
    nbl->ci[nbl->nci].cj_ind_end    = nbl->ncj;
}

/* Make a new sci entry at index nbl->nsci */
static void new_sci_entry(nbnxn_pairlist_t *nbl,int sci,int shift,int flags,
                          nbnxn_list_work_t *work)
{
    if (nbl->nsci + 1 > nbl->sci_nalloc)
    {
        nb_realloc_sci(nbl,nbl->nsci+1);
    }
    nbl->sci[nbl->nsci].sci           = sci;
    nbl->sci[nbl->nsci].shift         = shift;
    nbl->sci[nbl->nsci].cj4_ind_start = nbl->ncj4;
    nbl->sci[nbl->nsci].cj4_ind_end   = nbl->ncj4;
}

/* Sort the simple j-list cj on exclusions.
 * Entries with exclusions will all be sorted to the beginning of the list.
 */
static void sort_cj_excl(nbnxn_cj_t *cj,int ncj,
                         nbnxn_list_work_t *work)
{
    int jnew,j;

    if (ncj > work->cj_nalloc)
    {
        work->cj_nalloc = over_alloc_large(ncj);
        srenew(work->cj,work->cj_nalloc);
    }

    /* Make a list of the j-cells involving exclusions */
    jnew = 0;
    for(j=0; j<ncj; j++)
    {
        if (cj[j].excl != NBNXN_INT_MASK_ALL)
        {
            work->cj[jnew++] = cj[j];
        }
    }
    /* Check if there are exclusions at all or not just the first entry */
    if (!((jnew == 0) ||
          (jnew == 1 && cj[0].excl != NBNXN_INT_MASK_ALL)))
    {
        for(j=0; j<ncj; j++)
        {
            if (cj[j].excl == NBNXN_INT_MASK_ALL)
            {
                work->cj[jnew++] = cj[j];
            }
        }
        for(j=0; j<ncj; j++)
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
        sort_cj_excl(nbl->cj+nbl->ci[nbl->nci].cj_ind_start,jlen,nbl->work);

        if (nbl->ci[nbl->nci].shift & NBNXN_CI_HALF_LJ(0))
        {
            nbl->work->ncj_hlj += jlen;
        }
        else if (!(nbl->ci[nbl->nci].shift & NBNXN_CI_DO_COUL(0)))
        {
            nbl->work->ncj_noq += jlen;
        }

        nbl->nci++;
    }
}

/* Split sci entry for load balancing on the GPU.
 * As we only now the current count on our own thread,
 * we will need to estimate the current total amount of i-entries.
 * As the lists get concatenated later, this estimate depends
 * both on nthread and our own thread index thread.
 */
static void split_sci_entry(nbnxn_pairlist_t *nbl,
                            int nsp_max_av,gmx_bool progBal,int nc_bal,
                            int thread,int nthread)
{
    int nsci_est;
    int nsp_max;
    int cj4_start,cj4_end,j4len,cj4;
    int sci;
    int nsp,nsp_sci,nsp_cj4,nsp_cj4_e,nsp_cj4_p;
    int p;

    /* Estimate the total numbers of ci's of the nblist combined
     * over all threads using the target number of ci's.
     */
    nsci_est = nc_bal*thread/nthread + nbl->nsci;
    if (progBal)
    {
        /* The first ci blocks should be larger, to avoid overhead.
         * The last ci blocks should be smaller, to improve load balancing.
         */
        nsp_max = max(1,
                      nsp_max_av*nc_bal*3/(2*(nsci_est - 1 + nc_bal)));
    }
    else
    {
        nsp_max = nsp_max_av;
    }

    cj4_start = nbl->sci[nbl->nsci-1].cj4_ind_start;
    cj4_end   = nbl->sci[nbl->nsci-1].cj4_ind_end;
    j4len = cj4_end - cj4_start;

    if (j4len > 1 && j4len*GPU_NSUBCELL*NBNXN_GPU_JGROUP_SIZE > nsp_max)
    {
        /* Remove the last ci entry and process the cj4's again */
        nbl->nsci -= 1;

        sci        = nbl->nsci;
        cj4        = cj4_start;
        nsp        = 0;
        nsp_sci    = 0;
        nsp_cj4_e  = 0;
        nsp_cj4    = 0;
        while (cj4 < cj4_end)
        {
            nsp_cj4_p = nsp_cj4;
            nsp_cj4   = 0;
            for(p=0; p<GPU_NSUBCELL*NBNXN_GPU_JGROUP_SIZE; p++)
            {
                nsp_cj4 += (nbl->cj4[cj4].imei[0].imask >> p) & 1;
            }
            nsp += nsp_cj4;

            if (nsp > nsp_max && nsp > nsp_cj4)
            {
                nbl->sci[sci].cj4_ind_end = cj4;
                sci++;
                nbl->nsci++;
                if (nbl->nsci+1 > nbl->sci_nalloc)
                {
                    nb_realloc_sci(nbl,nbl->nsci+1);
                }
                nbl->sci[sci].sci           = nbl->sci[nbl->nsci-1].sci;
                nbl->sci[sci].shift         = nbl->sci[nbl->nsci-1].shift;
                nbl->sci[sci].cj4_ind_start = cj4;
                nsp_sci   = nsp - nsp_cj4;
                nsp_cj4_e = nsp_cj4_p;
                nsp       = nsp_cj4;
            }

            cj4++;
        }

        /* Put the remaining cj4's in a new ci entry */
        nbl->sci[sci].cj4_ind_end = cj4_end;

        /* Possibly balance out the last two ci's
         * by moving the last cj4 of the second last ci.
         */
        if (nsp_sci - nsp_cj4_e >= nsp + nsp_cj4_e)
        {
            nbl->sci[sci-1].cj4_ind_end--;
            nbl->sci[sci].cj4_ind_start--;
        }

        sci++;
        nbl->nsci++;
    }
}

/* Clost this super/sub list i entry */
static void close_ci_entry_supersub(nbnxn_pairlist_t *nbl,
                                    int nsp_max_av,
                                    gmx_bool progBal,int nc_bal,
                                    int thread,int nthread)
{
    int j4len,tlen;
    int nb,b;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    j4len = nbl->sci[nbl->nsci].cj4_ind_end - nbl->sci[nbl->nsci].cj4_ind_start;
    if (j4len > 0)
    {
        /* We can only have complete blocks of 4 j-entries in a list,
         * so round the count up before closing.
         */
        nbl->ncj4         = ((nbl->work->cj_ind + 4-1) >> 2);
        nbl->work->cj_ind = nbl->ncj4*NBNXN_GPU_JGROUP_SIZE;

        nbl->nsci++;

        if (nsp_max_av > 0)
        {
            split_sci_entry(nbl,nsp_max_av,progBal,nc_bal,thread,nthread);
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

/* Sets a simple list i-cell bounding box, including PBC shift */
static void set_icell_bb_simple(const float *bb,int ci,
                                real shx,real shy,real shz,
                                float *bb_ci)
{
    int ia;

    ia = ci*NNBSBB_B;
    bb_ci[BBL_X] = bb[ia+BBL_X] + shx;
    bb_ci[BBL_Y] = bb[ia+BBL_Y] + shy;
    bb_ci[BBL_Z] = bb[ia+BBL_Z] + shz;
    bb_ci[BBU_X] = bb[ia+BBU_X] + shx;
    bb_ci[BBU_Y] = bb[ia+BBU_Y] + shy;
    bb_ci[BBU_Z] = bb[ia+BBU_Z] + shz;
}

/* Sets a super-cell and sub cell bounding boxes, including PBC shift */
static void set_icell_bb_supersub(const float *bb,int ci,
                                  real shx,real shy,real shz,
                                  float *bb_ci)
{
    int ia,m,i;

#ifdef NBNXN_BBXXXX
    ia = ci*(GPU_NSUBCELL>>STRIDE_8BB_2LOG)*NNBSBB_XXXX;
    for(m=0; m<(GPU_NSUBCELL>>STRIDE_8BB_2LOG)*NNBSBB_XXXX; m+=NNBSBB_XXXX)
    {
        for(i=0; i<STRIDE_8BB; i++)
        {
            bb_ci[m+0*STRIDE_8BB+i] = bb[ia+m+0*STRIDE_8BB+i] + shx;
            bb_ci[m+1*STRIDE_8BB+i] = bb[ia+m+1*STRIDE_8BB+i] + shy;
            bb_ci[m+2*STRIDE_8BB+i] = bb[ia+m+2*STRIDE_8BB+i] + shz;
            bb_ci[m+3*STRIDE_8BB+i] = bb[ia+m+3*STRIDE_8BB+i] + shx;
            bb_ci[m+4*STRIDE_8BB+i] = bb[ia+m+4*STRIDE_8BB+i] + shy;
            bb_ci[m+5*STRIDE_8BB+i] = bb[ia+m+5*STRIDE_8BB+i] + shz;
        }
    }
#else
    ia = ci*GPU_NSUBCELL*NNBSBB_B;
    for(i=0; i<GPU_NSUBCELL*NNBSBB_B; i+=NNBSBB_B)
    {
        bb_ci[BBL_X] = bb[ia+BBL_X] + shx;
        bb_ci[BBL_Y] = bb[ia+BBL_Y] + shy;
        bb_ci[BBL_Z] = bb[ia+BBL_Z] + shz;
        bb_ci[BBU_X] = bb[ia+BBU_X] + shx;
        bb_ci[BBU_Y] = bb[ia+BBU_Y] + shy;
        bb_ci[BBU_Z] = bb[ia+BBU_Z] + shz;
    }
#endif
}

/* Copies PBC shifted i-cell atom coordinates x,y,z to working array */
static void icell_set_x_simple(int ci,
                               real shx,real shy,real shz,
                               int na_c,
                               int stride,const real *x,
                               nbnxn_list_work_t *work)
{
    int  ia,i;

    ia = ci*NBNXN_CPU_CLUSTER_I_SIZE;

    for(i=0; i<NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        work->x_ci[i*STRIDE_XYZ+XX] = x[(ia+i)*stride+XX] + shx;
        work->x_ci[i*STRIDE_XYZ+YY] = x[(ia+i)*stride+YY] + shy;
        work->x_ci[i*STRIDE_XYZ+ZZ] = x[(ia+i)*stride+ZZ] + shz;
    }
}

/* Copies PBC shifted super-cell atom coordinates x,y,z to working array */
static void icell_set_x_supersub(int ci,
                                 real shx,real shy,real shz,
                                 int na_c,
                                 int stride,const real *x,
                                 nbnxn_list_work_t *work)
{
    int  ia,i;
    real *x_ci;

    x_ci = work->x_ci;

    ia = ci*GPU_NSUBCELL*na_c;
    for(i=0; i<GPU_NSUBCELL*na_c; i++)
    {
        x_ci[i*DIM + XX] = x[(ia+i)*stride + XX] + shx;
        x_ci[i*DIM + YY] = x[(ia+i)*stride + YY] + shy;
        x_ci[i*DIM + ZZ] = x[(ia+i)*stride + ZZ] + shz;
    }
}

#ifdef NBNXN_SEARCH_SSE
/* Copies PBC shifted super-cell packed atom coordinates to working array */
static void icell_set_x_supersub_sse8(int ci,
                                      real shx,real shy,real shz,
                                      int na_c,
                                      int stride,const real *x,
                                      nbnxn_list_work_t *work)
{
    int  si,io,ia,i,j;
    real *x_ci;

    x_ci = work->x_ci;

    for(si=0; si<GPU_NSUBCELL; si++)
    {
        for(i=0; i<na_c; i+=STRIDE_8BB)
        {
            io = si*na_c + i;
            ia = ci*GPU_NSUBCELL*na_c + io;
            for(j=0; j<STRIDE_8BB; j++)
            {
                x_ci[io*DIM + j + XX*STRIDE_8BB] = x[(ia+j)*stride+XX] + shx;
                x_ci[io*DIM + j + YY*STRIDE_8BB] = x[(ia+j)*stride+YY] + shy;
                x_ci[io*DIM + j + ZZ*STRIDE_8BB] = x[(ia+j)*stride+ZZ] + shz;
            }
        }
    }
}
#endif

static real nbnxn_rlist_inc_nonloc_fac = 0.6;

/* Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 */
real nbnxn_get_rlist_effective_inc(int cluster_size,real atom_density)
{
    return ((0.5 + nbnxn_rlist_inc_nonloc_fac)*sqr(((cluster_size) - 1.0)/(cluster_size))*pow((cluster_size)/(atom_density),1.0/3.0));
}

/* Estimates the interaction volume^2 for non-local interactions */
static real nonlocal_vol2(const gmx_domdec_zones_t *zones,rvec ls,real r)
{
    int  z,d;
    real cl,ca,za;
    real vold_est;
    real vol2_est_tot;

    vol2_est_tot = 0;

    /* Here we simply add up the volumes of 1, 2 or 3 1D decomposition
     * not home interaction volume^2. As these volumes are not additive,
     * this is an overestimate, but it would only be significant in the limit
     * of small cells, where we anyhow need to split the lists into
     * as small parts as possible.
     */

    for(z=0; z<zones->n; z++)
    {
        if (zones->shift[z][XX] + zones->shift[z][YY] + zones->shift[z][ZZ] == 1)
        {
            cl = 0;
            ca = 1;
            za = 1;
            for(d=0; d<DIM; d++)
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
static int get_nsubpair_max(const nbnxn_search_t nbs,
                            int iloc,
                            real rlist,
                            int min_ci_balanced)
{
    const nbnxn_grid_t *grid;
    rvec ls;
    real xy_diag2,r_eff_sup,vol_est,nsp_est,nsp_est_nl;
    int  nsubpair_max;

    grid = &nbs->grid[0];

    ls[XX] = (grid->c1[XX] - grid->c0[XX])/(grid->ncx*GPU_NSUBCELL_X);
    ls[YY] = (grid->c1[YY] - grid->c0[YY])/(grid->ncy*GPU_NSUBCELL_Y);
    ls[ZZ] = (grid->c1[ZZ] - grid->c0[ZZ])*grid->ncx*grid->ncy/(grid->nc*GPU_NSUBCELL_Z);

    /* The average squared length of the diagonal of a sub cell */
    xy_diag2 = ls[XX]*ls[XX] + ls[YY]*ls[YY] + ls[ZZ]*ls[ZZ];

    /* The formulas below are a heuristic estimate of the average nsj per si*/
    r_eff_sup = rlist + nbnxn_rlist_inc_nonloc_fac*sqr((grid->na_c - 1.0)/grid->na_c)*sqrt(xy_diag2/3);

    if (!nbs->DomDec || nbs->zones->n == 1)
    {
        nsp_est_nl = 0;
    }
    else
    {
        nsp_est_nl =
            sqr(grid->atom_density/grid->na_c)*
            nonlocal_vol2(nbs->zones,ls,r_eff_sup);
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
        vol_est += 0.5*4.0/3.0*M_PI*pow(r_eff_sup,3);

        nsp_est = grid->nsubc_tot*vol_est*grid->atom_density/grid->na_c;

        /* Subtract the non-local pair count */
        nsp_est -= nsp_est_nl;

        if (debug)
        {
            fprintf(debug,"nsp_est local %5.1f non-local %5.1f\n",
                    nsp_est,nsp_est_nl);
        }
    }
    else
    {
        nsp_est = nsp_est_nl;
    }

    if (min_ci_balanced <= 0 || grid->nc >= min_ci_balanced || grid->nc == 0)
    {
        /* We don't need to worry */
        nsubpair_max = -1;
    }
    else
    {
        /* Thus the (average) maximum j-list size should be as follows */
        nsubpair_max = max(1,(int)(nsp_est/min_ci_balanced+0.5));

        /* Since the target value is a maximum (this avoid high outliers,
         * which lead to load imbalance), not average, we get more lists
         * than we ask for (to compensate we need to add GPU_NSUBCELL*4/4).
         * But more importantly, the optimal GPU performance moves
         * to lower number of block for very small blocks.
         * To compensate we add the maximum pair count per cj4.
         */
        nsubpair_max += GPU_NSUBCELL*NBNXN_CPU_CLUSTER_I_SIZE;
    }

    if (debug)
    {
        fprintf(debug,"nbl nsp estimate %.1f, nsubpair_max %d\n",
                nsp_est,nsubpair_max);
    }

    return nsubpair_max;
}

/* Debug list print function */
static void print_nblist_ci_cj(FILE *fp,const nbnxn_pairlist_t *nbl)
{
    int i,j;

    for(i=0; i<nbl->nci; i++)
    {
        fprintf(fp,"ci %4d  shift %2d  ncj %3d\n",
                nbl->ci[i].ci,nbl->ci[i].shift,
                nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start);

        for(j=nbl->ci[i].cj_ind_start; j<nbl->ci[i].cj_ind_end; j++)
        {
            fprintf(fp,"  cj %5d  imask %x\n",
                    nbl->cj[j].cj,
                    nbl->cj[j].excl);
        }
    }
}

/* Debug list print function */
static void print_nblist_sci_cj(FILE *fp,const nbnxn_pairlist_t *nbl)
{
    int i,j4,j;

    for(i=0; i<nbl->nsci; i++)
    {
        fprintf(fp,"ci %4d  shift %2d  ncj4 %2d\n",
                nbl->sci[i].sci,nbl->sci[i].shift,
                nbl->sci[i].cj4_ind_end - nbl->sci[i].cj4_ind_start);

        for(j4=nbl->sci[i].cj4_ind_start; j4<nbl->sci[i].cj4_ind_end; j4++)
        {
            for(j=0; j<4; j++)
            {
                fprintf(fp,"  sj %5d  imask %x\n",
                        nbl->cj4[j4].cj[j],
                        nbl->cj4[j4].imei[0].imask);
            }
        }
    }
}

/* Combine pair lists *nbl generated on multiple threads nblc */
static void combine_nblists(int nnbl,nbnxn_pairlist_t **nbl,
                            nbnxn_pairlist_t *nblc)
{
    int nsci,ncj4,nexcl;
    int n,i;

    if (nblc->bSimple)
    {
        gmx_incons("combine_nblists does not support simple lists");
    }

    nsci  = nblc->nsci;
    ncj4  = nblc->ncj4;
    nexcl = nblc->nexcl;
    for(i=0; i<nnbl; i++)
    {
        nsci  += nbl[i]->nsci;
        ncj4  += nbl[i]->ncj4;
        nexcl += nbl[i]->nexcl;
    }

    if (nsci > nblc->sci_nalloc)
    {
        nb_realloc_sci(nblc,nsci);
    }
    if (ncj4 > nblc->cj4_nalloc)
    {
        nblc->cj4_nalloc = over_alloc_small(ncj4);
        nbnxn_realloc_void((void **)&nblc->cj4,
                           nblc->ncj4*sizeof(*nblc->cj4),
                           nblc->cj4_nalloc*sizeof(*nblc->cj4),
                           nblc->alloc,nblc->free);
    }
    if (nexcl > nblc->excl_nalloc)
    {
        nblc->excl_nalloc = over_alloc_small(nexcl);
        nbnxn_realloc_void((void **)&nblc->excl,
                           nblc->nexcl*sizeof(*nblc->excl),
                           nblc->excl_nalloc*sizeof(*nblc->excl),
                           nblc->alloc,nblc->free);
    }

    /* Each thread should copy its own data to the combined arrays,
     * as otherwise data will go back and forth between different caches.
     */
#pragma omp parallel for num_threads(gmx_omp_nthreads_get(emntPairsearch)) schedule(static)
    for(n=0; n<nnbl; n++)
    {
        int sci_offset;
        int cj4_offset;
        int ci_offset;
        int excl_offset;
        int i,j4;
        const nbnxn_pairlist_t *nbli;

        /* Determine the offset in the combined data for our thread */
        sci_offset  = nblc->nsci;
        cj4_offset  = nblc->ncj4;
        ci_offset   = nblc->nci_tot;
        excl_offset = nblc->nexcl;

        for(i=0; i<n; i++)
        {
            sci_offset  += nbl[i]->nsci;
            cj4_offset  += nbl[i]->ncj4;
            ci_offset   += nbl[i]->nci_tot;
            excl_offset += nbl[i]->nexcl;
        }

        nbli = nbl[n];

        for(i=0; i<nbli->nsci; i++)
        {
            nblc->sci[sci_offset+i]                = nbli->sci[i];
            nblc->sci[sci_offset+i].cj4_ind_start += cj4_offset;
            nblc->sci[sci_offset+i].cj4_ind_end   += cj4_offset;
        }

        for(j4=0; j4<nbli->ncj4; j4++)
        {
            nblc->cj4[cj4_offset+j4] = nbli->cj4[j4];
            nblc->cj4[cj4_offset+j4].imei[0].excl_ind += excl_offset;
            nblc->cj4[cj4_offset+j4].imei[1].excl_ind += excl_offset;
        }

        for(j4=0; j4<nbli->nexcl; j4++)
        {
            nblc->excl[excl_offset+j4] = nbli->excl[j4];
        }
    }

    for(n=0; n<nnbl; n++)
    {
        nblc->nsci    += nbl[n]->nsci;
        nblc->ncj4    += nbl[n]->ncj4;
        nblc->nci_tot += nbl[n]->nci_tot;
        nblc->nexcl   += nbl[n]->nexcl;
    }
}

/* Returns the next ci to be processes by our thread */
static gmx_bool next_ci(const nbnxn_grid_t *grid,
                        int conv,
                        int nth,int ci_block,
                        int *ci_x,int *ci_y,
                        int *ci_b,int *ci)
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
                                        real rlist,
                                        gmx_bool simple)
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
    real bbx,bby;
    real rbb2;

    bbx = 0.5*(gridi->sx + gridj->sx);
    bby = 0.5*(gridi->sy + gridj->sy);
    if (!simple)
    {
        bbx /= GPU_NSUBCELL_X;
        bby /= GPU_NSUBCELL_Y;
    }

    rbb2 = sqr(max(0,rlist - 0.5*sqrt(bbx*bbx + bby*bby)));

#ifndef GMX_DOUBLE
    return rbb2;
#else
    return (float)((1+GMX_FLOAT_EPS)*rbb2);
#endif
}

static int get_ci_block_size(const nbnxn_grid_t *gridi,
                             gmx_bool bDomDec, int nth,
                             gmx_bool *bFBufferFlag)
{
    const int ci_block_enum = 5;
    const int ci_block_denom = 11;
    const int ci_block_min_atoms = 16;
    int ci_block;

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
    if (!bDomDec || ci_block*3*nth > gridi->nc)
    {
        ci_block = (gridi->nc + nth - 1)/nth;
        /* With non-interleaved blocks it makes sense to flag which
         * part of the force output thread buffer we access.
         * We use bit flags, so we have to check if it fits.
         */
        *bFBufferFlag = (nth > 1 && nth <= sizeof(unsigned int)*8);
    }
    else
    {
        *bFBufferFlag = FALSE;
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
                                     int min_ci_balanced,
                                     int th,int nth,
                                     nbnxn_pairlist_t *nbl)
{
    int  na_cj_2log;
    matrix box;
    real rl2;
    float rbb2;
    int  d;
    int  ci_b,ci,ci_x,ci_y,ci_xy,cj;
    ivec shp;
    int  tx,ty,tz;
    int  shift;
    gmx_bool bMakeList;
    real shx,shy,shz;
    int  conv_i,cell0_i;
    const float *bb_i,*bbcz_i,*bbcz_j;
    const int *flags_i;
    real bx0,bx1,by0,by1,bz0,bz1;
    real bz1_frac;
    real d2cx,d2z,d2z_cx,d2z_cy,d2zx,d2zxy,d2xy;
    int  cxf,cxl,cyf,cyf_x,cyl;
    int  cx,cy;
    int  c0,c1,cs,cf,cl;
    int  ndistc;
    int  ncpcheck;
    int  gridj_flag_shift=0,cj_offset=0;
    unsigned *gridj_flag=NULL;
    int  ncj_old_i,ncj_old_j;

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
        init_grid_flags(&work->gridi_flags,gridi);
        init_grid_flags(&work->gridj_flags,gridj);

        /* To flag j-blocks for gridj, we need to convert j-clusters to flag blocks */
        gridj_flag_shift = 0;
        while ((nbl->na_cj<<gridj_flag_shift) < NBNXN_CELLBLOCK_SIZE*nbl->na_ci)
        {
            gridj_flag_shift++;
        }
        /* We will subtract the cell offset, which is not a multiple of the block size */
        cj_offset = ci_to_cj(get_2log(nbl->na_cj),gridj->cell0);

        gridj_flag = work->gridj_flags.flag;
    }

    copy_mat(nbs->box,box);

    rl2 = nbl->rlist*nbl->rlist;

    rbb2 = boundingbox_only_distance2(gridi,gridj,nbl->rlist,nbl->bSimple);

    if (debug)
    {
        fprintf(debug,"nbl bounding box only distance %f\n",sqrt(rbb2));
    }

    /* Set the shift range */
    for(d=0; d<DIM; d++)
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
        bb_i    = gridi->bb;
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
        fprintf(debug,"nbl nc_i %d col.av. %.1f ci_block %d\n",
                gridi->nc,gridi->nc/(double)(gridi->ncx*gridi->ncy),ci_block);
    }

    ndistc = 0;
    ncpcheck = 0;

    /* Initially ci_b and ci to 1 before where we want them to start,
     * as they will both be incremented in next_ci.
     */
    ci_b = -1;
    ci   = th*ci_block - 1;
    ci_x = 0;
    ci_y = 0;
    while (next_ci(gridi,conv_i,nth,ci_block,&ci_x,&ci_y,&ci_b,&ci))
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
                bx1 = bb_i[ci*NNBSBB_B+NNBSBB_C+XX];
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
        for (tz=-shp[ZZ]; tz<=shp[ZZ]; tz++)
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

            for (ty=-shp[YY]; ty<=shp[YY]; ty++)
            {
                shy = ty*box[YY][YY] + tz*box[ZZ][YY];

                if (nbl->bSimple)
                {
                    by0 = bb_i[ci*NNBSBB_B         +YY] + shy;
                    by1 = bb_i[ci*NNBSBB_B+NNBSBB_C+YY] + shy;
                }
                else
                {
                    by0 = gridi->c0[YY] + (ci_y  )*gridi->sy + shy;
                    by1 = gridi->c0[YY] + (ci_y+1)*gridi->sy + shy;
                }

                get_cell_range(by0,by1,
                               gridj->ncy,gridj->c0[YY],gridj->sy,gridj->inv_sy,
                               d2z_cx,rl2,
                               &cyf,&cyl);

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

                for (tx=-shp[XX]; tx<=shp[XX]; tx++)
                {
                    shift = XYZ2IS(tx,ty,tz);

#ifdef NBNXN_SHIFT_BACKWARD
                    if (gridi == gridj && shift > CENTRAL)
                    {
                        continue;
                    }
#endif

                    shx = tx*box[XX][XX] + ty*box[YY][XX] + tz*box[ZZ][XX];

                    if (nbl->bSimple)
                    {
                        bx0 = bb_i[ci*NNBSBB_B         +XX] + shx;
                        bx1 = bb_i[ci*NNBSBB_B+NNBSBB_C+XX] + shx;
                    }
                    else
                    {
                        bx0 = gridi->c0[XX] + (ci_x  )*gridi->sx + shx;
                        bx1 = gridi->c0[XX] + (ci_x+1)*gridi->sx + shx;
                    }

                    get_cell_range(bx0,bx1,
                                   gridj->ncx,gridj->c0[XX],gridj->sx,gridj->inv_sx,
                                   d2z_cy,rl2,
                                   &cxf,&cxl);

                    if (cxf > cxl)
                    {
                        continue;
                    }

                    if (nbl->bSimple)
                    {
                        new_ci_entry(nbl,cell0_i+ci,shift,flags_i[ci],
                                     nbl->work);
                    }
                    else
                    {
                        new_sci_entry(nbl,cell0_i+ci,shift,flags_i[ci],
                                      nbl->work);
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
                        set_icell_bb_simple(bb_i,ci,shx,shy,shz,
                                            nbl->work->bb_ci);
                    }
                    else
                    {
                        set_icell_bb_supersub(bb_i,ci,shx,shy,shz,
                                              nbl->work->bb_ci);
                    }

                    nbs->icell_set_x(cell0_i+ci,shx,shy,shz,
                                     gridi->na_c,nbat->xstride,nbat->x,
                                     nbl->work);

                    for(cx=cxf; cx<=cxl; cx++)
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

                        for(cy=cyf_x; cy<=cyl; cy++)
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
                                while(cf > c0 &&
                                      (bbcz_j[cf*NNBSBB_D+1] >= bz0 ||
                                       d2xy + sqr(bbcz_j[cf*NNBSBB_D+1] - bz0) < rl2))
                                {
                                    cf--;
                                }

                                /* Find the highest cell that can possibly
                                 * be within range.
                                 */
                                cl = cs;
                                while(cl < c1-1 &&
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
                                    for(k=c0; k<c1; k++)
                                    {
                                        if (box_dist2(bx0,bx1,by0,by1,bz0,bz1,
                                                      bb+k*NNBSBB_B) < rl2 &&
                                            k < cf)
                                        {
                                            cf = k;
                                        }
                                        if (box_dist2(bx0,bx1,by0,by1,bz0,bz1,
                                                      bb+k*NNBSBB_B) < rl2 &&
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
                                    cf = max(cf,ci);
#else
                                    if (shift == CENTRAL)
                                    {
                                        cf = max(cf,ci);
                                    }
#endif
                                }

                                if (cf <= cl)
                                {
                                    /* For f buffer flags with simple lists */
                                    ncj_old_j = nbl->ncj;

                                    switch (nb_kernel_type)
                                    {
                                    case nbk4x4_PlainC:
                                        check_subcell_list_space_simple(nbl,cl-cf+1);

                                        make_cluster_list_simple(gridj,
                                                                 nbl,ci,cf,cl,
                                                                 (gridi == gridj && shift == CENTRAL),
                                                                 nbat->x,
                                                                 rl2,rbb2,
                                                                 &ndistc);
                                        break;
#ifdef NBNXN_SEARCH_SSE
                                    case nbk4xN_X86_SIMD128:
                                        check_subcell_list_space_simple(nbl,ci_to_cj(na_cj_2log,cl-cf)+2);
                                        make_cluster_list_x86_simd128(gridj,
                                                                      nbl,ci,cf,cl,
                                                                      (gridi == gridj && shift == CENTRAL),
                                                                      nbat->x,
                                                                      rl2,rbb2,
                                                                      &ndistc);
                                        break;
#ifdef GMX_X86_AVX_256
                                    case nbk4xN_X86_SIMD256:
                                        check_subcell_list_space_simple(nbl,ci_to_cj(na_cj_2log,cl-cf)+2);
                                        make_cluster_list_x86_simd256(gridj,
                                                                      nbl,ci,cf,cl,
                                                                      (gridi == gridj && shift == CENTRAL),
                                                                      nbat->x,
                                                                      rl2,rbb2,
                                                                      &ndistc);
                                        break;
#endif
#endif
                                    case nbk8x8x8_PlainC:
                                    case nbk8x8x8_CUDA:
                                        check_subcell_list_space_supersub(nbl,cl-cf+1);
                                        for(cj=cf; cj<=cl; cj++)
                                        {
                                            make_cluster_list_supersub(nbs,gridi,gridj,
                                                                       nbl,ci,cj,
                                                                       (gridi == gridj && shift == CENTRAL && ci == cj),
                                                                       nbat->xstride,nbat->x,
                                                                       rl2,rbb2,
                                                                       &ndistc);
                                        }
                                        break;
                                    }
                                    ncpcheck += cl - cf + 1;

                                    if (bFBufferFlag && nbl->ncj > ncj_old_j)
                                    {
                                        int cbf,cbl,cb;

                                        cbf = (nbl->cj[ncj_old_j].cj - cj_offset) >> gridj_flag_shift;
                                        cbl = (nbl->cj[nbl->ncj-1].cj - cj_offset) >> gridj_flag_shift;
                                        for(cb=cbf; cb<=cbl; cb++)
                                        {
                                            gridj_flag[cb] = 1U<<th;
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
                    }
                    else
                    {
                        set_sci_top_excls(nbs,
                                          nbl,
                                          shift == CENTRAL && gridi == gridj,
                                          gridj->na_c_2log,
                                          &(nbl->sci[nbl->nsci]),
                                          excl);
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
                                                progBal,min_ci_balanced,
                                                th,nth);
                    }
                }
            }
        }

        if (bFBufferFlag && nbl->ncj > ncj_old_i)
        {
            work->gridi_flags.flag[ci>>NBNXN_CELLBLOCK_SIZE_2LOG] = 1U<<th;
        }
    }

    work->ndistc = ndistc;

    nbs_cycle_stop(&work->cc[enbsCCsearch]);

    if (debug)
    {
        fprintf(debug,"number of distance checks %d\n",ndistc);
        fprintf(debug,"ncpcheck %s %d\n",gridi==gridj ? "local" : "non-local",
                ncpcheck);

        if (nbl->bSimple)
        {
            print_nblist_statistics_simple(debug,nbl,nbs,rlist);
        }
        else
        {
            print_nblist_statistics_supersub(debug,nbl,nbs,rlist);
        }

    }
}

static void reduce_cellblock_flags(const nbnxn_search_t nbs,
                                   int nnbl,
                                   const nbnxn_grid_t *gridi,
                                   const nbnxn_grid_t *gridj)
{
    int nbl,cb;
    const unsigned *flag;

    if (gridi->cellblock_flags.bUse)
    {
        for(nbl=0; nbl<nnbl; nbl++)
        {
            flag = nbs->work[nbl].gridi_flags.flag;
            
            for(cb=0; cb<gridi->cellblock_flags.ncb; cb++)
            {
                gridi->cellblock_flags.flag[cb] |= flag[cb];
            }
        }
    }
    if (gridj->cellblock_flags.bUse)
    {
        for(nbl=0; nbl<nnbl; nbl++)
        {
            flag = nbs->work[nbl].gridj_flags.flag;

            for(cb=0; cb<gridj->cellblock_flags.ncb; cb++)
            {
                gridj->cellblock_flags.flag[cb] |= flag[cb];
            }
        }
    }
}

static void print_reduction_cost(const nbnxn_grid_t *grids,int ngrid,int nnbl)
{
    int g,c0,c,cb,nbl;
    const nbnxn_grid_t *grid;

    for(g=0; g<ngrid; g++)
    {
        grid = &grids[g];

        c0 = 0;
        if (grid->cellblock_flags.bUse)
        {
            c  = 0;
            for(cb=0; cb<grid->cellblock_flags.ncb; cb++)
            {
                for(nbl=0; nbl<nnbl; nbl++)
                {
                    if (grid->cellblock_flags.flag[cb] == 1)
                    {
                        c0++;
                    }
                    else if (grid->cellblock_flags.flag[cb] & (1U<<nbl))
                    {
                        c++;
                    }
                }
            }
        }
        else
        {
            c = nnbl*grid->cellblock_flags.ncb;
        }
        fprintf(debug,"nbnxn reduction buffers, grid %d: %d flag %d only buf. 0: %4.2f av. reduction: %4.2f\n",
                g,nnbl,grid->cellblock_flags.bUse,
                c0/(double)(grid->cellblock_flags.ncb),
                c/(double)(grid->cellblock_flags.ncb));
    }
}

/* Make a local or non-local pair-list, depending on iloc */
void nbnxn_make_pairlist(const nbnxn_search_t nbs,
                         const nbnxn_atomdata_t *nbat,
                         const t_blocka *excl,
                         real rlist,
                         int min_ci_balanced,
                         nbnxn_pairlist_set_t *nbl_list,
                         int iloc,
                         int nb_kernel_type,
                         t_nrnb *nrnb)
{
    nbnxn_grid_t *gridi,*gridj;
    int nzi,zi,zj0,zj1,zj;
    int nsubpair_max;
    int th;
    int nnbl;
    nbnxn_pairlist_t **nbl;
    int ci_block;
    gmx_bool CombineNBLists,bFBufferFlag;
    int np_tot,np_noq,np_hlj,nap;

    nnbl            = nbl_list->nnbl;
    nbl             = nbl_list->nbl;
    CombineNBLists  = nbl_list->bCombined;

    if (debug)
    {
        fprintf(debug,"ns making %d nblists\n", nnbl);
    }

    if (nbl_list->bSimple)
    {
        switch (nb_kernel_type)
        {
#ifdef NBNXN_SEARCH_SSE
        case nbk4xN_X86_SIMD128:
            nbs->icell_set_x = icell_set_x_x86_simd128;
            break;
#ifdef GMX_X86_AVX_256
        case nbk4xN_X86_SIMD256:
            nbs->icell_set_x = icell_set_x_x86_simd256;
            break;
#endif
#endif
        default:
            nbs->icell_set_x = icell_set_x_simple;
            break;
        }
    }
    else
    {
#ifdef NBNXN_SEARCH_SSE
        nbs->icell_set_x = icell_set_x_supersub_sse8;
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
        nsubpair_max = get_nsubpair_max(nbs,iloc,rlist,min_ci_balanced);
    }
    else
    {
        nsubpair_max = 0;
    }

    /* Clear all pair-lists */
    for(th=0; th<nnbl; th++)
    {
        clear_pairlist(nbl[th]);
    }

    for(zi=0; zi<nzi; zi++)
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
        for(zj=zj0; zj<zj1; zj++)
        {
            gridj = &nbs->grid[zj];

            if (debug)
            {
                fprintf(debug,"ns search grid %d vs %d\n",zi,zj);
            }

            nbs_cycle_start(&nbs->cc[enbsCCsearch]);

            if (nbl[0]->bSimple && !gridi->bSimple)
            {
                /* Hybrid list, determine blocking later */
                ci_block = 0;
                bFBufferFlag = FALSE;
            }
            else
            {
                ci_block = get_ci_block_size(gridi,nbs->DomDec,nnbl,
                                             &bFBufferFlag);
                if (CombineNBLists)
                {
                    bFBufferFlag = FALSE;
                }
            }
            if (debug != NULL)
            {
                fprintf(debug,"grid %d %d F buffer flags %d\n",
                        zi,zj,bFBufferFlag);
            }

#pragma omp parallel for num_threads(nnbl) schedule(static)
            for(th=0; th<nnbl; th++)
            {
                if (CombineNBLists && th > 0)
                {
                    clear_pairlist(nbl[th]);
                }

                /* Divide the i super cell equally over the nblists */
                nbnxn_make_pairlist_part(nbs,gridi,gridj,
                                         &nbs->work[th],nbat,excl,
                                         rlist,
                                         nb_kernel_type,
                                         ci_block,
                                         bFBufferFlag,
                                         nsubpair_max,
                                         (LOCAL_I(iloc) || nbs->zones->n <= 2),
                                         min_ci_balanced,
                                         th,nnbl,
                                         nbl[th]);
            }
            nbs_cycle_stop(&nbs->cc[enbsCCsearch]);

            np_tot = 0;
            np_noq = 0;
            np_hlj = 0;
            for(th=0; th<nnbl; th++)
            {
                inc_nrnb(nrnb,eNR_NBNXN_DIST2,nbs->work[th].ndistc);

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
            nap = nbl[0]->na_ci*nbl[0]->na_cj;
            nbl_list->natpair_ljq = (np_tot - np_noq)*nap - np_hlj*nap/2;
            nbl_list->natpair_lj  = np_noq*nap;
            nbl_list->natpair_q   = np_hlj*nap/2;

            if (CombineNBLists && nnbl > 1)
            {
                nbs_cycle_start(&nbs->cc[enbsCCcombine]);

                combine_nblists(nnbl-1,nbl+1,nbl[0]);

                nbs_cycle_stop(&nbs->cc[enbsCCcombine]);
            }

            if (bFBufferFlag)
            {
                reduce_cellblock_flags(nbs,nnbl,gridi,gridj);
            }
            else
            {
                gridi->cellblock_flags.bUse = FALSE;
                gridj->cellblock_flags.bUse = FALSE;
            }
        }
    }

    /*
    print_supersub_nsp("nsubpair",nbl[0],iloc);
    */

    /* Special performance logging stuff (env.var. GMX_NBNXN_CYCLE) */
    if (LOCAL_I(iloc))
    {
        nbs->search_count++;
    }
    if (nbs->print_cycles &&
        (!nbs->DomDec || (nbs->DomDec && !LOCAL_I(iloc))) &&
        nbs->search_count % 100 == 0)
    {
        nbs_cycle_print(stderr,nbs);
    }

    if (debug && (CombineNBLists && nnbl > 1))
    {
        if (nbl[0]->bSimple)
        {
            print_nblist_statistics_simple(debug,nbl[0],nbs,rlist);
        }
        else
        {
            print_nblist_statistics_supersub(debug,nbl[0],nbs,rlist);
        }
    }

    if (gmx_debug_at)
    {
        if (nbl[0]->bSimple)
        {
            print_nblist_ci_cj(debug,nbl[0]);
        }
        else
        {
            print_nblist_sci_cj(debug,nbl[0]);
        }

        print_reduction_cost(nbs->grid,nbs->ngrid,nnbl);
    }
}
