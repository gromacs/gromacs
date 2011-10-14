/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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

#ifdef GMX_THREAD_SHM_FDECOMP
#include <pthread.h> 
#endif

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "pbc.h"
#include "nsbox.h"
#include "gmx_cyclecounter.h"
#include "gmxfio.h"

#ifdef GMX_OPENMP
#include <omp.h>
#else
#include "no_omp.h"
#endif

#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
#include "gmx_sse2_single.h"

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef GMX_SSE4_1
#include <smmintrin.h>
#endif
#endif


/* Neighbor search box upper and lower bound in x,y,z.
 * Store this in 4 iso 3 reals for SSE.
 */
#define NNBSBB_C  4
#define NNBSBB_B  (2*NNBSBB_C)
/* Neighbor search box upper and lower bound in z only. */
#define NNBSBB_D  2

#define SIMD_WIDTH       4
#define SIMD_WIDTH_2LOG  2

/* Interaction mask for 4vs4 SIMD interactions: all and only above diagonal */
#define SIMD_INT_MASK_ALL   0xffff
#define SIMD_INT_MASK_DIAG  0x8ce


#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
#define GMX_NBS_BBXXXX

#define NNBSBB_XXXX      (NNBSBB_D*DIM*SIMD_WIDTH)
#endif

#define NSBOX_SHIFT_BACKWARD


/* This define is a lazy way to avoid interdependence of the grid
 * and searching data structures.
 */
#define NBL_NAPC_MAX (NSUBCELL*16)
/* Working data for the actual i-supercell during neighbor searching */

typedef struct gmx_nbl_work {
    gmx_cache_protect_t cp0;

    real *bb_ci;       /* The bounding boxes, pbc shifted, for each subcell */
    real *x_ci;        /* The coordinates, pbc shifted, for each atom       */
    int  sj_ind;       /* The current sj_ind index for the current list     */
    int  sj4_init;     /* The first unitialized sj4 block                   */

    float *d2;         /* Bounding box distance work array                  */

#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
    __m128 ix_SSE0,iy_SSE0,iz_SSE0;
    __m128 ix_SSE1,iy_SSE1,iz_SSE1;
    __m128 ix_SSE2,iy_SSE2,iz_SSE2;
    __m128 ix_SSE3,iy_SSE3,iz_SSE3;
#endif

    gmx_nbl_cj_t *cj;
    int  cj_nalloc;

    gmx_cache_protect_t cp1;
} gmx_nbl_work_t;

typedef void
gmx_icell_set_x_t(int ci,
                  real shx,real shy,real shz,
                  int naps,
                  int stride,const real *x,
                  gmx_nbl_work_t *work);

static gmx_icell_set_x_t icell_set_x_simple;
static gmx_icell_set_x_t icell_set_x_simple_xxxx;
static gmx_icell_set_x_t icell_set_x_supersub;
static gmx_icell_set_x_t icell_set_x_supersub_sse8;

typedef gmx_bool
gmx_subcell_in_range_t(int naps,
                       int si,const real *x_or_bb_i,
                       int csj,int stride,const real *x_or_bb_j,
                       real rl2);

static gmx_subcell_in_range_t subc_in_range_x;
static gmx_subcell_in_range_t subc_in_range_sse8;

typedef struct {
    int          count;
    gmx_cycles_t c;
    gmx_cycles_t start;
} gmx_nbs_cc_t;

enum { enbsCCgrid, enbsCCsearch, enbsCCcombine, enbsCCreducef, enbsCCnr };

typedef struct {
    rvec c0;    /* The lower corner of the (local) grid         */
    rvec c1;    /* The upper corner of the (local) grid         */
    real atom_density;

    gmx_bool simple;
    int  na_tx;      /* Number of atoms per tixel            */
    int  na_stx;     /* Number of atoms per super-tixel      */
    int  na_tx_2log; /* 2log of na_tx                        */

    int  ncx;        /* Number of (super-)cells along x      */
    int  ncy;        /* Number of (super-)cells along y      */
    int  nc;         /* Total number of (super-)cells/tixels */

    real sx;         /* x-size of a (super-)cell             */
    real sy;         /* y-size of a (super-)cell             */
    real inv_sx;     /* 1/sx                                 */
    real inv_sy;     /* 1/sy                                 */

    int  cell0; /* The index in the nbs->cell array corresponding to cell 0 */

    int  *cxy_na;
    int  *cxy_ind;
    int  cxy_nalloc;

    int  *nsubc; /* The number of sub cells for each super cell */
    real *bbcz;  /* Bounding boxes in z for the super cells     */
    real *bb;    /* 3D bounding boxes for the sub cells         */
    int  *flags; /* Flag for the super cells                    */
    int  nc_nalloc;

    gmx_bool *use_simple;
    real *bbcz_simple;
    real *bb_simple;
    int  *flags_simple;
    int  nc_nalloc_simple;

    int  nsubc_tot;
} gmx_nbs_grid_t;

typedef struct {
    gmx_cache_protect_t cp0;

    int *cxy_na;
    int cxy_na_nalloc;

    int  *sort_work;
    int  sort_work_nalloc;

    float *bb_tmp;

    gmx_nbs_cc_t cc[enbsCCnr];

    gmx_cache_protect_t cp1;
} gmx_nbs_work_t;

typedef struct gmx_nbsearch {
    int  ePBC;
    matrix box;

    gmx_bool DomDec;
    ivec dd_dim;
    gmx_domdec_zones_t *zones;

    int  na_tx_supersub;      /* Number of atoms per tixel for non-simple */
    int  na_tx_supersub_2log; /* 2log of na_tx_supersub                   */

    int  ngrid;
    gmx_nbs_grid_t *grid;
    int  *cell; /* Actual allocated cell arry for all the grids */
    int  cell_nalloc;
    int  *a;    /* The atom index on the grid, the inverse of cell */
    int  a_nalloc;

    int  natoms_local;    /* The local atoms run from 0 to natoms_local */
    int  natoms_nonlocal; /* The non-local atoms run from natoms_local
                           * to natoms_nonlocal */

    gmx_bool print_cycles;
    gmx_nbs_cc_t cc[enbsCCnr];

    gmx_icell_set_x_t *icell_set_x;

    gmx_subcell_in_range_t *subc_dc;

    int  nthread_max;
    gmx_nbs_work_t *work;
} gmx_nbsearch_t_t;

static void nbs_cycle_clear(gmx_nbs_cc_t *cc)
{
    int i;

    for(i=0; i<enbsCCnr; i++)
    {
        cc[i].count = 0;
        cc[i].c     = 0;
    }
}

static void nbs_cycle_start(gmx_nbs_cc_t *cc)
{
    cc->start = gmx_cycles_read();
}

static void nbs_cycle_stop(gmx_nbs_cc_t *cc)
{
    cc->c += gmx_cycles_read() - cc->start;
    cc->count++;
}

static double Mcyc_av(const gmx_nbs_cc_t *cc)
{
    return (double)cc->c*1e-6/cc->count;
}

static void nbs_cycle_print(FILE *fp,const gmx_nbsearch_t nbs)
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
            fprintf(fp," comb %4.1f",
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

static void gmx_nbs_grid_init(gmx_nbs_grid_t * grid)
{
    grid->cxy_na      = NULL;
    grid->cxy_ind     = NULL;
    grid->cxy_nalloc  = 0;
    grid->bb          = NULL;
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
        gmx_fatal(FARGS,"nsbox naps (%d) is not a power of 2",n);
    }

    return log2;
}

void gmx_nbsearch_init(gmx_nbsearch_t * nbs_ptr,
                       ivec *n_dd_cells,
                       gmx_domdec_zones_t *zones,
                       gmx_bool simple,
                       int natoms_subcell,
                       int nthread_max)
{
    gmx_nbsearch_t nbs;
    int d,np,g,t;

    snew(nbs,1);
    *nbs_ptr = nbs;

    nbs->DomDec = (n_dd_cells != NULL);

    clear_ivec(nbs->dd_dim);
    if (!nbs->DomDec)
    {
        nbs->ngrid = 1;
    }
    else
    {
        nbs->zones = zones;

        np = 1;
        for(d=0; d<DIM; d++)
        {
            if ((*n_dd_cells)[d] > 1)
            {
                nbs->dd_dim[d] = 1;
                np *= 3;
            }
        }
        nbs->ngrid = (np + 1)/2;
    }

    if (!simple)
    {
        nbs->na_tx_supersub = natoms_subcell;

        if (nbs->na_tx_supersub > NBL_NAPC_MAX)
        {
            gmx_fatal(FARGS,
                      "napc (%d) is larger than the maximum allowed value (%d)",
                      nbs->na_tx_supersub,NBL_NAPC_MAX);
        }
        nbs->na_tx_supersub_2log = get_2log(nbs->na_tx_supersub);
    }

    snew(nbs->grid,nbs->ngrid);
    for(g=0; g<nbs->ngrid; g++)
    {
        gmx_nbs_grid_init(&nbs->grid[g]);
    }
    nbs->cell        = NULL;
    nbs->cell_nalloc = 0;
    nbs->a           = NULL;
    nbs->a_nalloc    = 0;

    /* nbs->subc_dc is only used with super/sub setup */
#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
    nbs->subc_dc = subc_in_range_sse8;
#else
    if (getenv("GMX_NSBOX_BB") != NULL)
    {
        /* Use only bounding box sub cell pair distances,
         * fast, but produces slightly more sub cell pairs.
         */
        nbs->subc_dc = NULL;
    }
    else
    {
        nbs->subc_dc = subc_in_range_x;
    }
#endif

    nbs->nthread_max = nthread_max;

    /* Initialize the work data structures for each thread */
    snew(nbs->work,nbs->nthread_max);
    for(t=0; t<nbs->nthread_max; t++)
    {
        nbs->work[t].cxy_na           = NULL;
        nbs->work[t].cxy_na_nalloc    = 0;
        nbs->work[t].sort_work        = NULL;
        nbs->work[t].sort_work_nalloc = 0;

        snew_aligned(nbs->work[t].bb_tmp,SIMD_WIDTH*NNBSBB_D,16);
    }

    /* Initialize detailed nbsearch cycle counting */
    nbs->print_cycles = (getenv("GMX_NBS_CYCLE") != 0);
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

static int set_grid_size_xy(const gmx_nbsearch_t nbs,
                            gmx_nbs_grid_t *grid,
                            int n,rvec corner0,rvec corner1,
                            real atom_density)
{
    rvec size;
    real adens,tlen,tlen_x,tlen_y,nc_max;
    int  t;

    rvec_sub(corner1,corner0,size);

    if (n > grid->na_stx)
    {
        /* target cell length */
        if (grid->simple)
        {
            /* Approximately cubic super cells */
            tlen   = pow(grid->na_stx/atom_density,1.0/3.0);
            tlen_x = tlen;
            tlen_y = tlen;
        }
        else
        {
            /* Approximately cubic sub cells */
            tlen   = pow(grid->na_tx/atom_density,1.0/3.0);
            tlen_x = tlen*NSUBCELL_X;
            tlen_y = tlen*NSUBCELL_Y;
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
    nc_max = n/grid->na_stx + grid->ncx*grid->ncy;
    if (nc_max > grid->nc_nalloc)
    {
        int bb_nalloc;

        grid->nc_nalloc = over_alloc_large(nc_max);
        srenew(grid->nsubc,grid->nc_nalloc);
        srenew(grid->bbcz,grid->nc_nalloc*NNBSBB_D);
#ifdef GMX_NBS_BBXXXX
        if (NSUBCELL % SIMD_WIDTH != 0)
        {
            gmx_incons("NSUBCELL is not a multiple of SIMD_WITH");
        }
        bb_nalloc = grid->nc_nalloc*NSUBCELL/SIMD_WIDTH*NNBSBB_XXXX;
#else  
        bb_nalloc = grid->nc_nalloc*NSUBCELL*NNBSBB_B;
#endif
        sfree_aligned(grid->bb);
        /* This snew also zeros the contents, this avoid possible
         * floating exceptions in SSE with the unused bb elements.
         */
        snew_aligned(grid->bb,bb_nalloc,16);

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

#ifdef DEBUG_NSBOX_GRIDDING
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

static void calc_bounding_box(int na,int stride,const real *x,real *bb)
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
    bb[0] = xl;
    bb[4] = xh;
    bb[1] = yl;
    bb[5] = yh;
    bb[2] = zl;
    bb[6] = zh;
}

/* Coordinate order xxxx, bb order xyz */
static void calc_bounding_box_x_xxxx(int na,const real *x,real *bb)
{
    int  j;
    real xl,xh,yl,yh,zl,zh;

    xl = x[0];
    xh = x[0];
    yl = x[4];
    yh = x[4];
    zl = x[8];
    zh = x[8];
    for(j=1; j<na; j++)
    {
        xl = min(xl,x[j+0]);
        xh = max(xh,x[j+0]);
        yl = min(yl,x[j+4]);
        yh = max(yh,x[j+4]);
        zl = min(zl,x[j+8]);
        zh = max(zh,x[j+8]);
    }
    bb[0] = xl;
    bb[4] = xh;
    bb[1] = yl;
    bb[5] = yh;
    bb[2] = zl;
    bb[6] = zh;
}

/* Coordinate order xyz, bb order xxxx */
static void calc_bounding_box_xxxx(int na,int stride,const real *x,real *bb)
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
    bb[ 0] = xl;
    bb[ 4] = yl;
    bb[ 8] = zl;
    bb[12] = xh;
    bb[16] = yh;
    bb[20] = zh;
}

#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )

static void calc_bounding_box_sse(int na,const real *x,float *bb)
{
    __m128 bb_0_SSE,bb_1_SSE;
    __m128 x_SSE;

    int  i;

    bb_0_SSE = _mm_load_ps(x);
    bb_1_SSE = bb_0_SSE;

    for(i=1; i<na; i++)
    {
        x_SSE    = _mm_load_ps(x+i*4);
        bb_0_SSE = _mm_min_ps(bb_0_SSE,x_SSE);
        bb_1_SSE = _mm_max_ps(bb_1_SSE,x_SSE);
    }

    _mm_store_ps(bb  ,bb_0_SSE);
    _mm_store_ps(bb+4,bb_1_SSE);
}

static void calc_bounding_box_xxxx_sse(int na,const real *x,
                                       float *bb_tmp,
                                       real *bb)
{
    calc_bounding_box_sse(na,x,bb_tmp);

    bb[ 0] = bb_tmp[0];
    bb[ 4] = bb_tmp[1];
    bb[ 8] = bb_tmp[2];
    bb[12] = bb_tmp[4];
    bb[16] = bb_tmp[5];
    bb[20] = bb_tmp[6];
}

#endif

static void print_bbsizes_simple(FILE *fp,
                                 const gmx_nbsearch_t nbs,
                                 const gmx_nbs_grid_t *grid)
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

static void print_bbsizes_supersub(FILE *fp,
                                   const gmx_nbsearch_t nbs,
                                   const gmx_nbs_grid_t *grid)
{
    int  ns,c,s;
    dvec ba;

    clear_dvec(ba);
    ns = 0;
    for(c=0; c<grid->nc; c++)
    {
#ifdef GMX_NBS_BBXXXX
        for(s=0; s<grid->nsubc[c]; s+=SIMD_WIDTH)
        {
            int cs_w,i,d;

            cs_w = (c*NSUBCELL + s)/SIMD_WIDTH;
            for(i=0; i<SIMD_WIDTH; i++)
            {
                for(d=0; d<DIM; d++)
                {
                    ba[d] +=
                        grid->bb[cs_w*NNBSBB_XXXX+(DIM+d)*SIMD_WIDTH+i] -
                        grid->bb[cs_w*NNBSBB_XXXX+     d *SIMD_WIDTH+i];
                }
            }
        }
#else
        for(s=0; s<grid->nsubc[c]; s++)
        {
            int cs,d;

            cs = c*NSUBCELL + s;
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
            nbs->box[XX][XX]/(grid->ncx*NSUBCELL_X),
            nbs->box[YY][YY]/(grid->ncy*NSUBCELL_Y),
            nbs->box[ZZ][ZZ]*grid->ncx*grid->ncy/(grid->nc*NSUBCELL_Z),
            ba[XX],ba[YY],ba[ZZ],
            ba[XX]*grid->ncx*NSUBCELL_X/nbs->box[XX][XX],
            ba[YY]*grid->ncy*NSUBCELL_Y/nbs->box[YY][YY],
            ba[ZZ]*grid->nc*NSUBCELL_Z/(grid->ncx*grid->ncy*nbs->box[ZZ][ZZ]));
}

static void copy_int_to_nbat_int(const int *a,int na,int na_round,
                                 const int *in,int fill,int *innb)
{
    int i,j;

    j = 0;
    for(i=0; i<na; i++)
    {
        innb[j++] = in[a[i]];
    }
    /* Complete the partially filled last cell with fill */
    for(; i<na_round; i++)
    {
        innb[j++] = fill;
    }
}

static void clear_nbat_real(int na,int stride,real *xnb)
{
    int i;

    for(i=0; i<na*stride; i++)
    {
        xnb[i] = 0;
    }
}

static void copy_rvec_to_nbat_real(const int *a,int na,int na_round,
                                   rvec *x,int nbatXFormat,real *xnb,
                                   int cx,int cy,int cz)
{
    int i,j,c;

/* We might need to place filler particles to fill ub the cell to na_round.
 * The coefficients (LJ and q) for such particles are zero.
 * But we might still get NaN as 0*NaN when distances are too small.
 * We hope that -107 nm is far away enough from to zero
 * to avoid accidental short distances to particles shifted down for pbc.
 */
#define NBAT_FAR_AWAY 107

    switch (nbatXFormat)
    {
    case nbatXYZ:
        j = 0;
        for(i=0; i<na; i++)
        {
            xnb[j++] = x[a[i]][XX];
            xnb[j++] = x[a[i]][YY];
            xnb[j++] = x[a[i]][ZZ];
        }
        /* Complete the partially filled last cell with copies of the last element.
         * This simplifies the bounding box calculation and avoid
         * numerical issues with atoms that are coincidentally close.
         */
        for(; i<na_round; i++)
        {
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cx);
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cy);
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cz + i);
        }
        break;
    case nbatXYZQ:
        j = 0;
        for(i=0; i<na; i++)
        {
            xnb[j++] = x[a[i]][XX];
            xnb[j++] = x[a[i]][YY];
            xnb[j++] = x[a[i]][ZZ];
            j++;
        }
        /* Complete the partially filled last cell with particles far apart */
        for(; i<na_round; i++)
        {
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cx);
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cy);
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cz + i);
            j++;
        }
        break;
    case nbatXXXX:
        j = 0;
        c = 0;
        for(i=0; i<na; i++)
        {
            xnb[j+0] = x[a[i]][XX];
            xnb[j+4] = x[a[i]][YY];
            xnb[j+8] = x[a[i]][ZZ];
            j++;
            c++;
            if (c == 4)
            {
                j += 8;
                c  = 0;
            }
        }
        /* Complete the partially filled last cell with particles far apart */
        for(; i<na_round; i++)
        {
            xnb[j+0] = -NBAT_FAR_AWAY*(1 + cx);
            xnb[j+4] = -NBAT_FAR_AWAY*(1 + cy);
            xnb[j+8] = -NBAT_FAR_AWAY*(1 + cz + i);
            j++;
            c++;
            if (c == 4)
            {
                j += 8;
                c  = 0;
            }
        }
        break;
    default:
        gmx_incons("Unsupported stride");
    }
}

void sort_on_lj(gmx_nb_atomdata_t *nbat,int naps,
                int a0,int a1,const int *atinfo,
                int *order,
                int *flags)
{
    int subc,s,a,ni,nj,a_lj_max,i,j;
    int sorti[NBL_NAPC_MAX/NSUBCELL];
    int sortj[NBL_NAPC_MAX/NSUBCELL];
    gmx_bool haveQ;

    subc = 0;
    for(s=a0; s<a1; s+=naps)
    {
        /* Sort this (sub-)cell on atoms with LH first, without last */
        ni = 0;
        nj = 0;
        haveQ = FALSE;
        a_lj_max = -1;
        for(a=s; a<min(s+naps,a1); a++)
        {
            haveQ = haveQ || GET_CGINFO_HAS_Q(atinfo[order[a]]);

            if (GET_CGINFO_HAS_LJ(atinfo[order[a]]))
            {
                sorti[ni++] = order[a];
                a_lj_max = a;
            }
            else
            {
                sortj[nj++] = order[a];
            }
        }

        *flags = 0;

        if (2*ni <= naps)
        {
            /* Only sort when strictly necessary, so we avoid summation
             * order precision loss as much as possible.
             */
            if (2*(a_lj_max - s) >= naps)
            {
                for(i=0; i<ni; i++)
                {
                    order[a0+i] = sorti[i];
                }
                for(j=0; j<nj; j++)
                {
                    order[a0+ni+j] = sortj[j];
                }
            }

            *flags |= NBL_CI_HALF_LJ(subc);
        }
        if (haveQ)
        {
            *flags |= NBL_CI_DO_COUL(subc);
        }
        subc++;
    }
}

void fill_cell(const gmx_nbsearch_t nbs,
               gmx_nbs_grid_t *grid,
               gmx_nb_atomdata_t *nbat,
               int a0,int a1,
               const int *atinfo,
               rvec *x,
               int sx,int sy, int sz)
{
    int  na,a;
    real *xnb;
    real *bb_ptr;

    na = a1 - a0;
    xnb = nbat->x + a0*nbat->xstride;

    if (na <= 0)
    {
        /* We only need to clear the dummy atom coordinates */
        clear_nbat_real(grid->na_tx,nbat->xstride,xnb);

        return;
    }

    if (grid->simple)
    {
        sort_on_lj(nbat,grid->na_tx,a0,a1,atinfo,nbs->a,
                   grid->flags+(a0>>grid->na_tx_2log)-grid->cell0);
    }

    /* Now we have sorted the atoms, set the cell indices */
    for(a=a0; a<a1; a++)
    {
        nbs->cell[nbs->a[a]] = a;
    }

    copy_rvec_to_nbat_real(nbs->a+a0,a1-a0,grid->na_tx,x,
                           nbat->XFormat,xnb,
                           sx,sy,sz);

    if (nbat->XFormat == nbatXXXX)
    {
        /* Store the bounding boxes as xyz.xyz. */
        bb_ptr =
            grid->bb+((a0-grid->cell0*grid->na_stx)>>grid->na_tx_2log)*NNBSBB_B;
        
        calc_bounding_box_x_xxxx(na,xnb,bb_ptr);
    }
#ifdef GMX_NBS_BBXXXX
    else if (!grid->simple)
    {
        /* Store the bounding boxes in a format convenient
         * for SSE calculations: xxxxyyyyzzzz...
                             */
        bb_ptr =
            grid->bb +
            ((a0-grid->cell0*grid->na_stx)>>(grid->na_tx_2log+SIMD_WIDTH_2LOG))*NNBSBB_XXXX +
            (((a0-grid->cell0*grid->na_stx)>>grid->na_tx_2log) & (SIMD_WIDTH-1));
        
        /* There is something wrong with this
         * calc_bounding_box_xxxx_sse function or call.
         * the results are incorrect.
         if (nbat->XFormat == nbatXYZQ)
         {
         calc_bounding_box_xxxx_sse(na,xnb,
                                                            nbs->work->bb_tmp,
                                                            bb_ptr);
                             }
                             else
        */
        {
            calc_bounding_box_xxxx(na,nbat->xstride,xnb,
                                   bb_ptr);
        }
        if (gmx_debug_at)
        {
            fprintf(debug,"%2d %2d %2d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    sx,sy,sz,
                    bb_ptr[0],bb_ptr[12],
                    bb_ptr[4],bb_ptr[16],
                    bb_ptr[8],bb_ptr[20]);
        }
    }
#endif
    else
    {
        /* Store the bounding boxes as xyz.xyz. */
        bb_ptr = grid->bb+((a0-grid->cell0*grid->na_stx)>>grid->na_tx_2log)*NNBSBB_B;
        
        calc_bounding_box(na,nbat->xstride,xnb,
                          bb_ptr);
        
        if (gmx_debug_at)
        {
            int bbo;
            bbo = (a0 - grid->cell0*grid->na_stx)/grid->na_tx;
            fprintf(debug,"%2d %2d %2d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    sx,sy,sz,
                    (grid->bb+bbo*NNBSBB_B)[0],
                    (grid->bb+bbo*NNBSBB_B)[4],
                    (grid->bb+bbo*NNBSBB_B)[1],
                    (grid->bb+bbo*NNBSBB_B)[5],
                    (grid->bb+bbo*NNBSBB_B)[2],
                    (grid->bb+bbo*NNBSBB_B)[6]);
        }
    }
}

static void sort_columns_simple(const gmx_nbsearch_t nbs,
                                int dd_zone,
                                gmx_nbs_grid_t *grid,
                                int a0,int a1,
                                const int *atinfo,
                                rvec *x,
                                gmx_nb_atomdata_t *nbat,
                                int cxy_start,int cxy_end,
                                int *sort_work)
{
    int  cxy;
    int  cx,cy,cz,ncz,c;
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
        ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_stx;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,FALSE,
                   nbs->a+ash,na,x,
                   grid->c0[ZZ],
                   ncz*grid->na_stx*SORT_GRID_OVERSIZE/nbs->box[ZZ][ZZ],
                   ncz*grid->na_stx*SGSF,sort_work);

        /* This loop goes over the supercells and subcells along z at once */
        for(cz=0; cz<ncz; cz++)
        {
            c  = grid->cxy_ind[cxy] + cz ;

            ash_c = ash + cz*grid->na_stx;
            na_c  = min(grid->na_stx,na-(ash_c-ash));
            
            fill_cell(nbs,grid,nbat,
                      ash_c,ash_c+na_c,atinfo,x,
                      grid->na_stx*cx + (dd_zone >> 2),
                      grid->na_stx*cy + (dd_zone & 3),
                      grid->na_stx*cz);

            /* This copy to bbcz is not really necessary.
             * But it allows to use the same grid search code
             * for the simple and supersub cell setups.
             */
            grid->bbcz[c*NNBSBB_D  ] = grid->bb[c*NNBSBB_B+2];
            grid->bbcz[c*NNBSBB_D+1] = grid->bb[c*NNBSBB_B+6];
        }

        /* Set the unused atom indices to -1 */
        for(ind=na; ind<ncz*grid->na_stx; ind++)
        {
            nbs->a[ash+ind] = -1;
        }
    }
}

static void sort_columns_supersub(const gmx_nbsearch_t nbs,
                                  int dd_zone,
                                  gmx_nbs_grid_t *grid,
                                  int a0,int a1,
                                  const int *atinfo,
                                  rvec *x,
                                  gmx_nb_atomdata_t *nbat,
                                  int cxy_start,int cxy_end,
                                  int *sort_work)
{
    int  cxy;
    int  cx,cy,cz=-1,c=-1,ncz;
    int  na,ash,na_c,ind,a;
    int  subdiv_z,sub_z,na_z,ash_z;
    int  subdiv_y,sub_y,na_y,ash_y;
    int  subdiv_x,sub_x,na_x,ash_x;

    if (debug)
    {
        fprintf(debug,"cell0 %d sorting columns %d - %d, atoms %d - %d\n",
                grid->cell0,cxy_start,cxy_end,a0,a1);
    }

    subdiv_x = grid->na_tx;
    subdiv_y = NSUBCELL_X*subdiv_x;
    subdiv_z = NSUBCELL_Y*subdiv_y;

    /* Sort the atoms within each x,y column in 3 dimensions */
    for(cxy=cxy_start; cxy<cxy_end; cxy++)
    {
        cx = cxy/grid->ncy;
        cy = cxy - cx*grid->ncy;

        na  = grid->cxy_na[cxy];
        ncz = grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy];
        ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_stx;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,FALSE,
                   nbs->a+ash,na,x,
                   grid->c0[ZZ],
                   ncz*grid->na_stx*SORT_GRID_OVERSIZE/nbs->box[ZZ][ZZ],
                   ncz*grid->na_stx*SGSF,sort_work);

        /* This loop goes over the supercells and subcells along z at once */
        for(sub_z=0; sub_z<ncz*NSUBCELL_Z; sub_z++)
        {
            ash_z = ash + sub_z*subdiv_z;
            na_z  = min(subdiv_z,na-(ash_z-ash));

            /* We have already sorted on z */

            if (sub_z % NSUBCELL_Z == 0)
            {
                cz = sub_z/NSUBCELL_Z;
                c  = grid->cxy_ind[cxy] + cz ;

                /* The number of atoms in this supercell */
                na_c = min(grid->na_stx,na-(ash_z-ash));

                grid->nsubc[c] = min(NSUBCELL,(na_c+grid->na_tx-1)/grid->na_tx);

                /* Store the z-boundaries of the super cell */
                grid->bbcz[c*NNBSBB_D  ] = x[nbs->a[ash_z]][ZZ];
                grid->bbcz[c*NNBSBB_D+1] = x[nbs->a[ash_z+na_c-1]][ZZ];
            }

#if NSUBCELL_Y > 1
            sort_atoms(YY,(sub_z & 1),
                       nbs->a+ash_z,na_z,x,
                       grid->c0[YY]+cy*grid->sy,grid->inv_sy,
                       subdiv_y*SGSF,sort_work);
#endif

            for(sub_y=0; sub_y<NSUBCELL_Y; sub_y++)
            {
                ash_y = ash_z + sub_y*subdiv_y;
                na_y  = min(subdiv_y,na-(ash_y-ash));

#if NSUBCELL_X > 1
                sort_atoms(XX,((cz*NSUBCELL_Y + sub_y) & 1),
                           nbs->a+ash_y,na_y,x,
                           grid->c0[XX]+cx*grid->sx,grid->inv_sx,
                           subdiv_x*SGSF,sort_work);
#endif

                for(sub_x=0; sub_x<NSUBCELL_X; sub_x++)
                {
                    ash_x = ash_y + sub_x*subdiv_x;
                    na_x  = min(subdiv_x,na-(ash_x-ash));

                    fill_cell(nbs,grid,nbat,
                              ash_x,ash_x+na_x,atinfo,x,
                              grid->na_tx*(cx*NSUBCELL_X+sub_x) + (dd_zone >> 2),
                              grid->na_tx*(cy*NSUBCELL_Y+sub_y) + (dd_zone & 3),
                              grid->na_tx*sub_z);
                }
            }
        }

        /* Set the unused atom indices to -1 */
        for(ind=na; ind<ncz*grid->na_stx; ind++)
        {
            nbs->a[ash+ind] = -1;
        }

        /*
        copy_rvec_to_nbat_real(axy,na,ncz*grid->na_stx,x,nbat->xstride,xnb);

        calc_bounding_box(ncz,grid->na_stx,nbat->xstride,xnb,
                          grid->bb+grid->cxy_ind[i]*NNBSBB);
        */
    }
}

static void calc_cell_indices(const gmx_nbsearch_t nbs,
                              int dd_zone,
                              gmx_nbs_grid_t *grid,
                              int a0,int a1,
                              const int *atinfo,
                              rvec *x,
                              int *move,
                              gmx_nb_atomdata_t *nbat)
{
    int  n0,n1,i;
    int  cx,cy,cxy,ncz_max,ncz;
    int  nthread,t;
    int  *cxy_na,cxy_na_i;

    nthread = omp_get_max_threads();

#pragma omp parallel  for private(cxy_na,n0,n1,i,cx,cy)
    for(t=0; t<nthread; t++)
    {
        cxy_na = nbs->work[t].cxy_na;

        /* We add one extra cell for particles which moved during DD */
        for(i=0; i<grid->ncx*grid->ncy+1; i++)
        {
            cxy_na[i] = 0;
        }

        n0 = a0 + (int)((t+0)*(a1 - a0))/nthread;
        n1 = a0 + (int)((t+1)*(a1 - a0))/nthread;
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
                nbs->cell[i] = cx*grid->ncy + cy;

#ifdef DEBUG_NSBOX_GRIDDING
                if (nbs->cell[i] < 0 || nbs->cell[i] >= grid->ncx*grid->ncy)
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
                nbs->cell[i] = grid->ncx*grid->ncy;
            }
            
            cxy_na[nbs->cell[i]]++;
        }
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
        for(t=1; t<nthread; t++)
        {
            cxy_na_i += nbs->work[t].cxy_na[i];
        }
        ncz = (cxy_na_i + grid->na_stx - 1)/grid->na_stx;
        grid->cxy_ind[i+1] = grid->cxy_ind[i] + ncz;
        /* Clear cxy_na, so we can reuse the array below */
        grid->cxy_na[i] = 0;
    }
    grid->nc = grid->cxy_ind[grid->ncx*grid->ncy] - grid->cxy_ind[0];

    nbat->natoms = (grid->cell0 + grid->nc)*grid->na_stx;

    if (debug)
    {
        fprintf(debug,"ns napc %d naps %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                grid->na_stx,grid->na_tx,grid->nc,
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
    if (ncz_max*grid->na_stx*SGSF > nbs->work[0].sort_work_nalloc)
    {
        for(t=0; t<nbs->nthread_max; t++)
        {
            nbs->work[t].sort_work_nalloc =
                over_alloc_large(ncz_max*grid->na_stx*SGSF);
            srenew(nbs->work[t].sort_work,nbs->work[t].sort_work_nalloc);
        }
    }

    /* Now we know the dimensions we can fill the grid.
     * This is the first, unsorted fill. We sort the columns after this.
     */
    for(i=a0; i<a1; i++)
    {
        /* At this point nbs->cell contains the local grid x,y indices */
        cxy = nbs->cell[i];
        nbs->a[(grid->cell0 + grid->cxy_ind[cxy])*grid->na_stx + grid->cxy_na[cxy]++] = i;
    }

    /* Set the cell indices for the moved particles */
    n0 = grid->nc*grid->na_stx;
    n1 = grid->nc*grid->na_stx+grid->cxy_na[grid->ncx*grid->ncy];
    for(i=n0; i<n1; i++)
    {
        nbs->cell[nbs->a[i]] = i;
    }

    /* Sort the super-cell columns along z into the sub-cells. */
    nthread = omp_get_max_threads();
#pragma omp parallel for schedule(static)
    for(t=0; t<nthread; t++)
    {
        if (grid->simple)
        {
            sort_columns_simple(nbs,dd_zone,grid,a0,a1,atinfo,x,nbat,
                                ((t+0)*grid->ncx*grid->ncy)/nthread,
                                ((t+1)*grid->ncx*grid->ncy)/nthread,
                                nbs->work[t].sort_work);
        }
        else
        {
            sort_columns_supersub(nbs,dd_zone,grid,a0,a1,atinfo,x,nbat,
                                  ((t+0)*grid->ncx*grid->ncy)/nthread,
                                  ((t+1)*grid->ncx*grid->ncy)/nthread,
                                  nbs->work[t].sort_work);
        }
    }

    if (!grid->simple)
    {
        grid->nsubc_tot = 0;
        for(i=0; i<grid->nc; i++)
        {
            grid->nsubc_tot += grid->nsubc[i];
        }
    }

    if (debug)
    {
        if (grid->simple)
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

static void nb_realloc_void(void **ptr,
                            int nbytes_copy,int nbytes_new,
                            gmx_nbat_alloc_t *ma,
                            gmx_nbat_free_t  *mf)
{
    void *ptr_new;

    ma(&ptr_new,nbytes_new);

    if (nbytes_new > 0 && ptr_new == NULL)
    {
        gmx_fatal(FARGS, "Allocation of %d bytes failed", nbytes_new);
    }

    if (nbytes_copy > 0)
    {
        if (nbytes_new < nbytes_copy)
        {
            gmx_incons("In nb_realloc_void: new size less than copy size");
        }
        memcpy(ptr_new,*ptr,nbytes_copy);
    }
    if (*ptr != NULL)
    {
        mf(*ptr);
    }
    *ptr = ptr_new;
}

/* NOTE: does not preserve the contents! */
static void nb_realloc_int(int **ptr,int n,
                           gmx_nbat_alloc_t *ma,
                           gmx_nbat_free_t  *mf)
{
    if (*ptr != NULL)
    {
        mf(*ptr);
    }
    ma((void **)ptr,n*sizeof(**ptr));
}

/* NOTE: does not preserve the contents! */
static void nb_realloc_real(real **ptr,int n,
                            gmx_nbat_alloc_t *ma,
                            gmx_nbat_free_t  *mf)
{
    if (*ptr != NULL)
    {
        mf(*ptr);
    }
    ma((void **)ptr,n*sizeof(**ptr));
}

static void gmx_nb_atomdata_realloc(gmx_nb_atomdata_t *nbat,int n)
{
    int t;

    nb_realloc_void((void **)&nbat->type,
                    nbat->natoms*sizeof(*nbat->type),
                    n*sizeof(*nbat->type),
                    nbat->alloc,nbat->free);
    nb_realloc_void((void **)&nbat->lj_comb,
                    nbat->natoms*2*sizeof(*nbat->lj_comb),
                    n*2*sizeof(*nbat->lj_comb),
                    nbat->alloc,nbat->free);
    if (nbat->XFormat != nbatXYZQ)
    {
        nb_realloc_void((void **)&nbat->q,
                        nbat->natoms*sizeof(*nbat->q),
                        n*sizeof(*nbat->q),
                        nbat->alloc,nbat->free);
    }
    if (nbat->nenergrp > 1)
    {
        nb_realloc_void((void **)&nbat->energrp,
                        nbat->natoms/nbat->naps*sizeof(*nbat->energrp),
                        n/nbat->naps*sizeof(*nbat->energrp),
                        nbat->alloc,nbat->free);
    }
    nb_realloc_void((void **)&nbat->x,
                    nbat->natoms*nbat->xstride*sizeof(*nbat->x),
                    n*nbat->xstride*sizeof(*nbat->x),
                    nbat->alloc,nbat->free);
    for(t=0; t<nbat->nout; t++)
    {
        nb_realloc_void((void **)&nbat->out[t].f,
                        nbat->natoms*nbat->xstride*sizeof(*nbat->out[t].f),
                        n*nbat->xstride*sizeof(*nbat->out[t].f),
                        nbat->alloc,nbat->free);
    }
    nbat->nalloc = n;
}

void gmx_nbsearch_put_on_grid(gmx_nbsearch_t nbs,
                              int ePBC,matrix box,
                              int dd_zone,
                              rvec corner0,rvec corner1,
                              int a0,int a1,
                              real atom_density,
                              const int *atinfo,
                              rvec *x,
                              int nmoved,int *move,
                              gmx_bool simple,
                              gmx_nb_atomdata_t *nbat)
{
    gmx_nbs_grid_t *grid;
    int n;
    int nc_max_grid,nc_max;

    grid = &nbs->grid[dd_zone];

    nbs_cycle_start(&nbs->cc[enbsCCgrid]);

    grid->simple = simple;

    if (grid->simple)
    {
        grid->na_tx      = SIMD_WIDTH;
        grid->na_stx     = SIMD_WIDTH;
        grid->na_tx_2log = SIMD_WIDTH_2LOG;
    }
    else
    {
        grid->na_tx      = nbs->na_tx_supersub;
        grid->na_stx     = NSUBCELL*nbs->na_tx_supersub;
        grid->na_tx_2log = nbs->na_tx_supersub_2log;
    }

    nbat->naps = grid->na_tx;

    if (dd_zone == 0)
    {
        grid->cell0 = 0;
    }
    else
    {
        grid->cell0 =
            (nbs->grid[dd_zone-1].cell0 + nbs->grid[dd_zone-1].nc)*
            nbs->grid[dd_zone-1].na_stx/grid->na_stx;
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
        /* We assume that gmx_nbsearch_put_on_grid is called first
         * for the local atoms (dd_zone=0).
         */
        nbs->natoms_nonlocal = a1 - nmoved;
    }
    else
    {
        nbs->natoms_nonlocal = max(nbs->natoms_nonlocal,a1);
    }

    nc_max_grid = set_grid_size_xy(nbs,grid,n-nmoved,corner0,corner1,
                                   nbs->grid[0].atom_density);

    nc_max = grid->cell0 + nc_max_grid;

    if (a1 > nbs->cell_nalloc)
    {
        nbs->cell_nalloc = over_alloc_large(a1);
        srenew(nbs->cell,nbs->cell_nalloc);
    }

    if (nc_max*grid->na_stx > nbs->a_nalloc)
    {
        nbs->a_nalloc = over_alloc_large(nc_max*grid->na_stx);
        srenew(nbs->a,nbs->a_nalloc);
    }

    if (nc_max*grid->na_stx > nbat->nalloc)
    {
        gmx_nb_atomdata_realloc(nbat,nc_max*grid->na_stx);
    }
    
    calc_cell_indices(nbs,dd_zone,grid,a0,a1,atinfo,x,move,nbat);

    if (dd_zone == 0)
    {
        nbat->natoms_local = nbat->natoms;
    }

    nbs_cycle_stop(&nbs->cc[enbsCCgrid]);
}

void gmx_nbsearch_put_on_grid_nonlocal(gmx_nbsearch_t nbs,
                                       const gmx_domdec_zones_t *zones,
                                       const int *atinfo,
                                       rvec *x,
                                       gmx_bool simple,
                                       gmx_nb_atomdata_t *nbat)
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

        gmx_nbsearch_put_on_grid(nbs,nbs->ePBC,NULL,
                                 zone,c0,c1,
                                 zones->cg_range[zone],
                                 zones->cg_range[zone+1],
                                 -1,
                                 atinfo,
                                 x,
                                 0,NULL,
                                 simple,
                                 nbat);
    }
}

void gmx_nbsearch_grid_simple(gmx_nbsearch_t nbs,
                              gmx_nb_atomdata_t *nbat)
{
    gmx_nbs_grid_t *grid;
    real *bbcz,*bb;
    int ncd,sc;

    grid = &nbs->grid[0];

    if (grid->simple)
    {
        gmx_incons("gmx_nbsearch_grid_simple called with a simple grid");
    }

    ncd = grid->na_stx/SIMD_WIDTH;

    if (grid->nc*ncd > grid->nc_nalloc_simple)
    {
        grid->nc_nalloc_simple = over_alloc_large(grid->nc*ncd);
        srenew(grid->use_simple,grid->nc_nalloc_simple);
        srenew(grid->bbcz_simple,grid->nc_nalloc_simple*NNBSBB_D);
        srenew(grid->bb_simple,grid->nc_nalloc_simple*NNBSBB_B);
        srenew(grid->flags_simple,grid->nc_nalloc_simple);
    }

    bbcz = grid->bbcz_simple;
    bb   = grid->bb_simple;

#pragma omp parallel for schedule(static)
    for(sc=0; sc<grid->nc; sc++)
    {
        int c,tx,na;

        for(c=0; c<ncd; c++)
        {
            tx = sc*ncd + c;

            na = SIMD_WIDTH;
            while(nbat->type[tx*SIMD_WIDTH+na-1] == nbat->ntype-1)
            {
                na--;
            }

            if (na > 0)
            {
                grid->use_simple[tx] = TRUE;

                if (nbat->XFormat == nbatXXXX)
                {
                    calc_bounding_box_x_xxxx(na,nbat->x+tx*SIMD_WIDTH*nbat->xstride,
                                             bb+tx*NNBSBB_B);
                }
                else
                {
                    calc_bounding_box(na,nbat->xstride,
                                      nbat->x+tx*SIMD_WIDTH*nbat->xstride,
                                      bb+tx*NNBSBB_B);
                }
                bbcz[tx*NNBSBB_D+0] = bb[tx*NNBSBB_B         +ZZ];
                bbcz[tx*NNBSBB_D+1] = bb[tx*NNBSBB_B+NNBSBB_C+ZZ];

                /* No interaction optimization yet here */
                grid->flags_simple[tx] = NBL_CI_DO_COUL(0);
            }
            else
            {
                grid->use_simple[tx] = FALSE;
            }
        }
    }
}

void gmx_nbsearch_get_ncells(gmx_nbsearch_t nbs,int *ncx,int *ncy)
{
    *ncx = nbs->grid[0].ncx;
    *ncy = nbs->grid[0].ncy;
}

void gmx_nbsearch_get_atomorder(gmx_nbsearch_t nbs,int **a,int *moved)
{
    const gmx_nbs_grid_t *grid;

    grid = &nbs->grid[0];

    /* Return the atom order for the home cell (index 0) */
    *a  = nbs->cell;

    *moved = grid->cxy_ind[grid->ncx*grid->ncy]*grid->na_stx;
}

void gmx_nbsearch_set_atomorder(gmx_nbsearch_t nbs)
{
    gmx_nbs_grid_t *grid;
    int ao,cx,cy,cxy,cz,j;

    /* Set the atom order for the home cell (index 0) */
    grid = &nbs->grid[0];

    ao = 0;
    for(cx=0; cx<grid->ncx; cx++)
    {
        for(cy=0; cy<grid->ncy; cy++)
        {
            cxy = cx*grid->ncy + cy;
            j   = grid->cxy_ind[cxy]*grid->na_stx;
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

static real box_dist2(real bx0,real bx1,real by0,real by1,real bz0,real bz1,
                      const real *bb)
{
    real d2;
    real dl,dh,dm,dm0;

    d2 = 0;

    dl  = bx0 - bb[4];
    dh  = bb[0] - bx1;
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = by0 - bb[5];
    dh  = bb[1] - by1;
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = bz0 - bb[6];
    dh  = bb[2] - bz1;
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    return d2;
}

static real subc_bb_dist2(int naps,
                          int si,const real *bb_i_ci,
                          int csj,const real *bb_j_all)
{
    const real *bb_i,*bb_j;
    real d2;
    real dl,dh,dm,dm0;

    bb_i = bb_i_ci  +  si*NNBSBB_B;
    bb_j = bb_j_all + csj*NNBSBB_B;

    d2 = 0;

    dl  = bb_i[0] - bb_j[4];
    dh  = bb_j[0] - bb_i[4];
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = bb_i[1] - bb_j[5];
    dh  = bb_j[1] - bb_i[5];
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    dl  = bb_i[2] - bb_j[6];
    dh  = bb_j[2] - bb_i[6];
    dm  = max(dl,dh);
    dm0 = max(dm,0);
    d2 += dm0*dm0;

    return d2;
}

#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )

static real subc_bb_dist2_sse(int naps,
                              int si,const real *bb_i_ci,
                              int csj,const real *bb_j_all)
{
    const real *bb_i,*bb_j;

    __m128 bb_i_SSE0,bb_i_SSE1;
    __m128 bb_j_SSE0,bb_j_SSE1;
    __m128 dl_SSE;
    __m128 dh_SSE;
    __m128 dm_SSE;
    __m128 dm0_SSE;
    __m128 d2_SSE;

    float d2_array[7],*d2_align;

    d2_align = (float *)(((size_t)(d2_array+3)) & (~((size_t)15)));

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
    d2_SSE    = _mm_mul_ps(dm0_SSE,dm0_SSE);

    _mm_store_ps(d2_align,d2_SSE);

    return d2_align[0] + d2_align[1] + d2_align[2];
}

static void subc_bb_dist2_sse_xxxx(const real *bb_j,
                                   int nsi,const real *bb_i,
                                   float *d2)
{
    int si;
    int shi;

    __m128 xj_l,yj_l,zj_l;
    __m128 xj_h,yj_h,zj_h;
    __m128 xi_l,yi_l,zi_l;
    __m128 xi_h,yi_h,zi_h;

    __m128 dx_0,dy_0,dz_0;
    __m128 dx_1,dy_1,dz_1;

    __m128 mx,my,mz;
    __m128 m0x,m0y,m0z;

    __m128 d2x,d2y,d2z;
    __m128 d2s,d2t;

    __m128 zero;

    zero = _mm_setzero_ps();

    xj_l = _mm_load1_ps(bb_j+0*SIMD_WIDTH);
    yj_l = _mm_load1_ps(bb_j+1*SIMD_WIDTH);
    zj_l = _mm_load1_ps(bb_j+2*SIMD_WIDTH);
    xj_h = _mm_load1_ps(bb_j+3*SIMD_WIDTH);
    yj_h = _mm_load1_ps(bb_j+4*SIMD_WIDTH);
    zj_h = _mm_load1_ps(bb_j+5*SIMD_WIDTH);

    for(si=0; si<nsi; si+=SIMD_WIDTH)
    {
        shi = si*NNBSBB_D*DIM;

        xi_l = _mm_load_ps(bb_i+shi+0*SIMD_WIDTH);
        yi_l = _mm_load_ps(bb_i+shi+1*SIMD_WIDTH);
        zi_l = _mm_load_ps(bb_i+shi+2*SIMD_WIDTH);
        xi_h = _mm_load_ps(bb_i+shi+3*SIMD_WIDTH);
        yi_h = _mm_load_ps(bb_i+shi+4*SIMD_WIDTH);
        zi_h = _mm_load_ps(bb_i+shi+5*SIMD_WIDTH);

        dx_0 = _mm_sub_ps(xi_l,xj_h);
        dy_0 = _mm_sub_ps(yi_l,yj_h);
        dz_0 = _mm_sub_ps(zi_l,zj_h);

        dx_1 = _mm_sub_ps(xj_l,xi_h);
        dy_1 = _mm_sub_ps(yj_l,yi_h);
        dz_1 = _mm_sub_ps(zj_l,zi_h);

        mx   = _mm_max_ps(dx_0,dx_1);
        my   = _mm_max_ps(dy_0,dy_1);
        mz   = _mm_max_ps(dz_0,dz_1);

        m0x  = _mm_max_ps(mx,zero);
        m0y  = _mm_max_ps(my,zero);
        m0z  = _mm_max_ps(mz,zero);

        d2x  = _mm_mul_ps(m0x,m0x);
        d2y  = _mm_mul_ps(m0y,m0y);
        d2z  = _mm_mul_ps(m0z,m0z);

        d2s  = _mm_add_ps(d2x,d2y);
        d2t  = _mm_add_ps(d2s,d2z);

        _mm_store_ps(d2+si,d2t);
    }
}

#ifdef GMX_SSE4_1
static real subc_bb_dist2_sse4_1(int naps,
                                 int si,const real *bb_i_ci,
                                 int csj,const real *bb_j_all)
{
    const real *bb_i,*bb_j;

    __m128 bb_i_SSE0,bb_i_SSE1;
    __m128 bb_j_SSE0,bb_j_SSE1;
    __m128 dl_SSE;
    __m128 dh_SSE;
    __m128 dm_SSE;
    __m128 dm0_SSE;
    __m128 d2_SSE;

    float d2;

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
    /* Dot product of components 0,1,2 */
    d2_SSE    = _mm_dp_ps(dm0_SSE,dm0_SSE,0x71);

    _mm_store_ss(&d2,d2_SSE);

    return d2;
}
#endif

#endif

static gmx_bool subc_in_range_x(int naps,
                                int si,const real *x_i,
                                int csj,int stride,const real *x_j,
                                real rl2)
{
    int  i,j,i0,j0;
    real d2;

    for(i=0; i<naps; i++)
    {
        i0 = (si*naps + i)*DIM;
        for(j=0; j<naps; j++)
        {
            j0 = (csj*naps + j)*stride;

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

static gmx_bool subc_in_range_sse8(int naps,
                                   int si,const real *x_i,
                                   int csj,int stride,const real *x_j,
                                   real rl2)
{
#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
    __m128 ix_SSE0,iy_SSE0,iz_SSE0;
    __m128 ix_SSE1,iy_SSE1,iz_SSE1;
    __m128 jx0_SSE,jy0_SSE,jz0_SSE;
    __m128 jx1_SSE,jy1_SSE,jz1_SSE;

    __m128     dx_SSE0,dy_SSE0,dz_SSE0;
    __m128     dx_SSE1,dy_SSE1,dz_SSE1;
    __m128     dx_SSE2,dy_SSE2,dz_SSE2;
    __m128     dx_SSE3,dy_SSE3,dz_SSE3;

    __m128     rsq_SSE0;
    __m128     rsq_SSE1;
    __m128     rsq_SSE2;
    __m128     rsq_SSE3;

    __m128     wco_SSE0;
    __m128     wco_SSE1;
    __m128     wco_SSE2;
    __m128     wco_SSE3;
    __m128     wco_any_SSE01,wco_any_SSE23,wco_any_SSE;
    
    __m128 rc2_SSE;

    float      wco_any_array[7],*wco_any_align;

    int naps_sse;
    int j0,j1;

    rc2_SSE   = _mm_set1_ps(rl2);

    wco_any_align = (float *)(((size_t)(wco_any_array+3)) & (~((size_t)15)));

    naps_sse = 8/4;
    ix_SSE0 = _mm_load_ps(x_i+(si*naps_sse*DIM+0)*4);
    iy_SSE0 = _mm_load_ps(x_i+(si*naps_sse*DIM+1)*4);
    iz_SSE0 = _mm_load_ps(x_i+(si*naps_sse*DIM+2)*4);
    ix_SSE1 = _mm_load_ps(x_i+(si*naps_sse*DIM+3)*4);
    iy_SSE1 = _mm_load_ps(x_i+(si*naps_sse*DIM+4)*4);
    iz_SSE1 = _mm_load_ps(x_i+(si*naps_sse*DIM+5)*4);

    j0 = csj*naps;
    j1 = j0 + naps - 1;
    while (j0 < j1)
    {
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

        _mm_store_ps(wco_any_align,wco_any_SSE);

        if (wco_any_align[0] != 0 ||
            wco_any_align[1] != 0 || 
            wco_any_align[2] != 0 ||
            wco_any_align[3] != 0)
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

static int nbl_sj(const gmx_nblist_t *nbl,int sj_ind)
{
    return nbl->sj4[sj_ind>>2].sj[sj_ind & 3];
}

static unsigned nbl_imask0(const gmx_nblist_t *nbl,int sj_ind)
{
    return nbl->sj4[sj_ind>>2].imei[0].imask;
}

static void check_excl_space(gmx_nblist_t *nbl,int extra)
{
    if (nbl->nexcl+extra > nbl->excl_nalloc)
    {
        nbl->excl_nalloc = over_alloc_small(nbl->nexcl+extra);
        nb_realloc_void((void **)&nbl->excl,
                        nbl->nexcl*sizeof(*nbl->excl),
                        nbl->excl_nalloc*sizeof(*nbl->excl),
                        nbl->alloc,nbl->free);
    }
}

static void check_subcell_list_space_simple(gmx_nblist_t *nbl,
                                            int ncell)
{
    int cj_max;

    cj_max = nbl->ncj + ncell;

    if (cj_max > nbl->cj_nalloc)
    {
        nbl->cj_nalloc = over_alloc_small(cj_max);
        nb_realloc_void((void **)&nbl->cj,
                        nbl->ncj*sizeof(*nbl->cj),
                        nbl->cj_nalloc*sizeof(*nbl->cj),
                        nbl->alloc,nbl->free);
    }
}

static void check_subcell_list_space_supersub(gmx_nblist_t *nbl,
                                              int nsupercell)
{
    int nsj4_max,j4,j,w,t;

#define NWARP       2
#define WARP_SIZE  32

    /* We can have maximally nsupercell*NSUBCELL sj lists */
    /* We can store 4 j-subcell - i-supercell pairs in one struct.
     * since we round down, we need one extra entry.
     */
    nsj4_max = ((nbl->work->sj_ind + nsupercell*NSUBCELL + 4-1) >> 2);

    if (nsj4_max > nbl->sj4_nalloc)
    {
        nbl->sj4_nalloc = over_alloc_small(nsj4_max);
        nb_realloc_void((void **)&nbl->sj4,
                        nbl->work->sj4_init*sizeof(*nbl->sj4),
                        nbl->sj4_nalloc*sizeof(*nbl->sj4),
                        nbl->alloc,nbl->free);
    }

    if (nsj4_max > nbl->work->sj4_init)
    {
        for(j4=nbl->work->sj4_init; j4<nsj4_max; j4++)
        {
            /* No i-subcells and no excl's in the list initially */
            for(w=0; w<NWARP; w++)
            {
                nbl->sj4[j4].imei[w].imask    = 0U;
                nbl->sj4[j4].imei[w].excl_ind = 0;

            }
        }
        nbl->work->sj4_init = nsj4_max;
    }
}

static void nblist_alloc_aligned(void **ptr,size_t nbytes)
{
    *ptr = save_calloc_aligned("ptr",__FILE__,__LINE__,nbytes,1,16,0);
}

static void nblist_free_aligned(void *ptr)
{
    sfree_aligned(ptr);
}

static void set_no_excls(gmx_nbl_excl_t *excl)
{
    int t;

    for(t=0; t<WARP_SIZE; t++)
    {
        /* Turn all interaction bits on */
        excl->pair[t] = 0xffffffffU;
    }
}

void gmx_nbl_list_init(gmx_nbl_lists_t *nbl_list,
                       gmx_bool simple, gmx_bool combined,
                       gmx_nbat_alloc_t *alloc,
                       gmx_nbat_free_t  *free)
{
    int i;

    nbl_list->simple    = simple;
    nbl_list->combined  = combined;

    nbl_list->nnbl = 1;
#ifdef GMX_OPENMP
        nbl_list->nnbl = omp_get_max_threads(); // FIXME: remove omp_get_max_threads
#endif
        snew(nbl_list->nbl,nbl_list->nnbl);
#pragma omp parallel for schedule(static)
        for(i=0; i<nbl_list->nnbl; i++)
        {
            /* Allocate the nblist data structure locally on each thread
             * to optimize memory access for NUMA architectures.
             */
            snew(nbl_list->nbl[i],1);
            gmx_nblist_init(nbl_list->nbl[i],
                            /* Only list 0 is used on the GPU */
                            (i==0) ? alloc : NULL,
                            (i==0) ? free  : NULL);
        }

}

void gmx_nblist_init(gmx_nblist_t *nbl,
                     gmx_nbat_alloc_t *alloc,
                     gmx_nbat_free_t  *free)
{
    if (alloc == NULL)
    {
        nbl->alloc = nblist_alloc_aligned;
    }
    else
    {
        nbl->alloc = alloc;
    }
    if (free == NULL)
    {
        nbl->free = nblist_free_aligned;
    }
    else
    {
        nbl->free = free;
    }

    nbl->napc        = 0;
    nbl->naps        = 0;
    nbl->nci         = 0;
    nbl->ci          = NULL;
    nbl->ci_nalloc   = 0;
    nbl->ncj         = 0;
    nbl->cj          = NULL;
    nbl->cj_nalloc   = 0;
    nbl->nsj4        = 0;
    /* We need one element extra in sj, so alloc initially with 1 */
    nbl->sj4_nalloc  = 0;
    nbl->sj4         = NULL;
    nbl->nsi         = 0;

    nbl->excl_nalloc = 0;
    nbl->nexcl       = 0;
    check_excl_space(nbl,1);
    nbl->nexcl       = 1;
    set_no_excls(&nbl->excl[0]);

    snew(nbl->work,1);
#ifdef GMX_NBS_BBXXXX
    snew_aligned(nbl->work->bb_ci,NSUBCELL/SIMD_WIDTH*NNBSBB_XXXX,16);
#else
    snew_aligned(nbl->work->bb_ci,NSUBCELL*NNBSBB_B,16);
#endif
    snew_aligned(nbl->work->x_ci,NBL_NAPC_MAX*DIM,16);
    snew_aligned(nbl->work->d2,NSUBCELL,16);
}

static void print_nblist_statistics_simple(FILE *fp,const gmx_nblist_t *nbl,
                                           const gmx_nbsearch_t nbs,real rl)
{
    const gmx_nbs_grid_t *grid;
    int cs[SHIFTS];
    int s,i,j;
    int npexcl;

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp,"nbl nci %d ncj %d\n",
            nbl->nci,nbl->ncj);
    fprintf(fp,"nbl napc %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->napc,rl,nbl->ncj,nbl->ncj/(double)grid->nc,
            nbl->ncj/(double)grid->nc*grid->na_stx,
            nbl->ncj/(double)grid->nc*grid->na_stx/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nc*grid->na_stx/det(nbs->box)));

    fprintf(fp,"nbl average j cell list length %.1f\n",
            0.25*nbl->ncj/(double)nbl->nci);

    for(s=0; s<SHIFTS; s++)
    {
        cs[s] = 0;
    }
    npexcl = 0;
    for(i=0; i<nbl->nci; i++)
    {
        cs[nbl->ci[i].shift & NBL_CI_SHIFT] +=
            nbl->ci[i].cj_ind_end - nbl->ci[i].cj_ind_start;

        j = nbl->ci[i].cj_ind_start;
        while (j < nbl->ci[i].cj_ind_end &&
               nbl->cj[j].excl != SIMD_INT_MASK_ALL)
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

static void print_nblist_statistics_supersub(FILE *fp,const gmx_nblist_t *nbl,
                                             const gmx_nbsearch_t nbs,real rl)
{
    const gmx_nbs_grid_t *grid;
    int i,j4,j,si,b;
    int c[NSUBCELL+1];

    /* This code only produces correct statistics with domain decomposition */
    grid = &nbs->grid[0];

    fprintf(fp,"nbl nci %d nsj4 %d nsi %d excl4 %d\n",
            nbl->nci,nbl->nsj4,nbl->nsi,nbl->nexcl);
    fprintf(fp,"nbl naps %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->naps,rl,nbl->nsi,nbl->nsi/(double)grid->nsubc_tot,
            nbl->nsi/(double)grid->nsubc_tot*grid->na_tx,
            nbl->nsi/(double)grid->nsubc_tot*grid->na_tx/(0.5*4.0/3.0*M_PI*rl*rl*rl*grid->nsubc_tot*grid->na_tx/det(nbs->box)));

    fprintf(fp,"nbl average j super cell list length %.1f\n",
            0.25*nbl->nsj4/(double)nbl->nci);
    fprintf(fp,"nbl average i sub cell list length %.1f\n",
            nbl->nsi/(0.25*nbl->nsj4));

    for(si=0; si<=NSUBCELL; si++)
    {
        c[si] = 0;
    }
    for(i=0; i<nbl->nci; i++)
    {
        for(j4=nbl->ci[i].sj4_ind_start; j4<nbl->ci[i].sj4_ind_end; j4++)
        {
            for(j=0; j<4; j++)
            {
                b = 0;
                for(si=0; si<NSUBCELL; si++)
                {
                    if (nbl->sj4[j4].imei[0].imask & (1U << (j*NSUBCELL + si)))
                    {
                        b++;
                    }
                }
                c[b]++;
            }
        }
    }
    for(b=0; b<=NSUBCELL; b++)
    {
        fprintf(fp,"nbl j-list #i-subcell %d %7d %4.1f\n",
                b,c[b],100.0*c[b]/(double)(nbl->nsj4*4));
    }
}

static void print_supersub_nsp(const char *fn,
                               const gmx_nblist_t *nbl,
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
        for(j4=nbl->ci[i].sj4_ind_start; j4<nbl->ci[i].sj4_ind_end; j4++)
        {
            for(p=0; p<4*NSUBCELL; p++)
            {
                nsp += (nbl->sj4[j4].imei[0].imask >> p) & 1;
            }
        }
        fprintf(fp,"%4d %3d %3d\n",
                i,
                nsp,
                nbl->ci[i].sj4_ind_end-nbl->ci[i].sj4_ind_start);
    }

    fclose(fp);
}

static void low_get_nbl_exclusions(gmx_nblist_t *nbl,int sj4,
                                   int warp,gmx_nbl_excl_t **excl)
{
    if (nbl->sj4[sj4].imei[warp].excl_ind == 0)
    {
        /* No exclusions set, make a new list entry */
        nbl->sj4[sj4].imei[warp].excl_ind = nbl->nexcl;
        nbl->nexcl++;
        *excl = &nbl->excl[nbl->sj4[sj4].imei[warp].excl_ind];
        set_no_excls(*excl);
    }
    else
    {
        /* We already have some exclusions, new ones can be added to the list */
        *excl = &nbl->excl[nbl->sj4[sj4].imei[warp].excl_ind];
    }
}

static void get_nbl_exclusions_1(gmx_nblist_t *nbl,int sj4,
                                 int warp,gmx_nbl_excl_t **excl)
{
    if (nbl->sj4[sj4].imei[warp].excl_ind == 0)
    {
        /* We need to make a new list entry, check if we have space */
        check_excl_space(nbl,1);
    }
    low_get_nbl_exclusions(nbl,sj4,warp,excl);
}

static void get_nbl_exclusions_2(gmx_nblist_t *nbl,int sj4,
                                 gmx_nbl_excl_t **excl_w0,
                                 gmx_nbl_excl_t **excl_w1)
{
    /* Check for space we might need */
    check_excl_space(nbl,2);
    
    low_get_nbl_exclusions(nbl,sj4,0,excl_w0);
    low_get_nbl_exclusions(nbl,sj4,1,excl_w1);
}

static void set_self_and_newton_excls(gmx_nblist_t *nbl,
                                      int sj4_ind,int sj_offset,
                                      int si)
{
    gmx_nbl_excl_t *excl[2];
    int  ei,ej,w;

    /* Here we only set the set self and double pair exclusions */

    get_nbl_exclusions_2(nbl,sj4_ind,&excl[0],&excl[1]);
    
    if (nbl->TwoWay)
    {
        /* Only minor != major bits set */
        for(ej=0; ej<nbl->naps; ej++)
        {
            w = (ej>>2);
            excl[w]->pair[(ej&(4-1))*nbl->naps+ej] &=
                ~(1U << (sj_offset*NSUBCELL+si));
        }
    }
    else
    {
        /* Only minor < major bits set */
        for(ej=0; ej<nbl->naps; ej++)
        {
            w = (ej>>2);
            for(ei=ej; ei<nbl->naps; ei++)
            {
                excl[w]->pair[(ej&(4-1))*nbl->naps+ei] &=
                    ~(1U << (sj_offset*NSUBCELL+si));
            }
        }
    }
}

static void make_subcell_list_simple(const gmx_nbs_grid_t *gridj,
                                     gmx_nblist_t *nbl,
                                     int ci,int cjf,int cjl,
                                     gmx_bool remove_sub_diag,
                                     const real *x_j,
                                     real rl2,real rbb2)
{
    const gmx_nbl_work_t *work;

    const real *bb_ci;
    const real *x_ci;

    gmx_bool   InRange;
    real       d2;
    int        cjf_gl,cjl_gl,cj;

    work = nbl->work;

    bb_ci = nbl->work->bb_ci;
    x_ci  = nbl->work->x_ci;

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
        d2 = subc_bb_dist2(4,0,bb_ci,cjf,gridj->bb);
        
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
            for(i=0; i<4 && !InRange; i++)
            {
                for(j=0; j<4; j++)
                {
                    InRange = InRange ||
                        (sqr(x_ci[i*DIM+0] - x_j[(cjf_gl*4+j)*DIM+0]) +
                         sqr(x_ci[i*DIM+1] - x_j[(cjf_gl*4+j)*DIM+1]) +
                         sqr(x_ci[i*DIM+2] - x_j[(cjf_gl*4+j)*DIM+2]) < rl2);
                }
            }
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
        d2 = subc_bb_dist2(4,0,bb_ci,cjl,gridj->bb);
        
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
            for(i=0; i<4 && !InRange; i++)
            {
                for(j=0; j<4; j++)
                {
                    InRange = InRange ||
                        (sqr(x_ci[i*DIM+0] - x_j[(cjl_gl*4+j)*DIM+0]) +
                         sqr(x_ci[i*DIM+1] - x_j[(cjl_gl*4+j)*DIM+1]) +
                         sqr(x_ci[i*DIM+2] - x_j[(cjl_gl*4+j)*DIM+2]) < rl2);
                }
            }
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
            nbl->cj[nbl->ncj].c    = gridj->cell0 + cj;
            nbl->cj[nbl->ncj].excl = (remove_sub_diag && ci == cj ?
                                      SIMD_INT_MASK_DIAG : SIMD_INT_MASK_ALL);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

static void make_subcell_list_simple_xxxx(const gmx_nbs_grid_t *gridj,
                                          gmx_nblist_t *nbl,
                                          int ci,int cjf,int cjl,
                                          gmx_bool remove_sub_diag,
                                          const real *x_j,
                                          real rl2,real rbb2)
{
#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )

    const gmx_nbl_work_t *work;

    const real *bb_ci;

    __m128 jx_SSE,jy_SSE,jz_SSE;

    __m128     dx_SSE0,dy_SSE0,dz_SSE0;
    __m128     dx_SSE1,dy_SSE1,dz_SSE1;
    __m128     dx_SSE2,dy_SSE2,dz_SSE2;
    __m128     dx_SSE3,dy_SSE3,dz_SSE3;

    __m128     rsq_SSE0;
    __m128     rsq_SSE1;
    __m128     rsq_SSE2;
    __m128     rsq_SSE3;

    __m128     wco_SSE0;
    __m128     wco_SSE1;
    __m128     wco_SSE2;
    __m128     wco_SSE3;
    __m128     wco_any_SSE01,wco_any_SSE23,wco_any_SSE;
    
    __m128 rc2_SSE;

    float      wco_any_array[7],*wco_any_align;

    gmx_bool   InRange;
    real       d2;
    int        cjf_gl,cjl_gl,cj;

    work = nbl->work;

    bb_ci = nbl->work->bb_ci;

    rc2_SSE   = _mm_set1_ps(rl2);

    wco_any_align = (float *)(((size_t)(wco_any_array+3)) & (~((size_t)15)));

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
#ifdef GMX_SSE4_1
        d2 = subc_bb_dist2_sse4_1(4,0,bb_ci,cjf,gridj->bb);
#else
        d2 = subc_bb_dist2_sse(4,0,bb_ci,cjf,gridj->bb);
#endif
        
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
            cjf_gl = gridj->cell0 + cjf;
            jx_SSE  = _mm_load_ps(x_j+cjf_gl*12+0);
            jy_SSE  = _mm_load_ps(x_j+cjf_gl*12+4);
            jz_SSE  = _mm_load_ps(x_j+cjf_gl*12+8);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(work->ix_SSE0,jx_SSE);
            dy_SSE0            = _mm_sub_ps(work->iy_SSE0,jy_SSE);
            dz_SSE0            = _mm_sub_ps(work->iz_SSE0,jz_SSE);
            dx_SSE1            = _mm_sub_ps(work->ix_SSE1,jx_SSE);
            dy_SSE1            = _mm_sub_ps(work->iy_SSE1,jy_SSE);
            dz_SSE1            = _mm_sub_ps(work->iz_SSE1,jz_SSE);
            dx_SSE2            = _mm_sub_ps(work->ix_SSE2,jx_SSE);
            dy_SSE2            = _mm_sub_ps(work->iy_SSE2,jy_SSE);
            dz_SSE2            = _mm_sub_ps(work->iz_SSE2,jz_SSE);
            dx_SSE3            = _mm_sub_ps(work->ix_SSE3,jx_SSE);
            dy_SSE3            = _mm_sub_ps(work->iy_SSE3,jy_SSE);
            dz_SSE3            = _mm_sub_ps(work->iz_SSE3,jz_SSE);
            
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
            
            _mm_store_ps(wco_any_align,wco_any_SSE);
            
            InRange = (wco_any_align[0] != 0 ||
                       wco_any_align[1] != 0 || 
                       wco_any_align[2] != 0 ||
                       wco_any_align[3] != 0);
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
#ifdef GMX_SSE4_1
        d2 = subc_bb_dist2_sse4_1(4,0,bb_ci,cjl,gridj->bb);
#else
        d2 = subc_bb_dist2_sse(4,0,bb_ci,cjl,gridj->bb);
#endif
        
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
            cjl_gl = gridj->cell0 + cjl;
            jx_SSE  = _mm_load_ps(x_j+cjl_gl*12+0);
            jy_SSE  = _mm_load_ps(x_j+cjl_gl*12+4);
            jz_SSE  = _mm_load_ps(x_j+cjl_gl*12+8);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(work->ix_SSE0,jx_SSE);
            dy_SSE0            = _mm_sub_ps(work->iy_SSE0,jy_SSE);
            dz_SSE0            = _mm_sub_ps(work->iz_SSE0,jz_SSE);
            dx_SSE1            = _mm_sub_ps(work->ix_SSE1,jx_SSE);
            dy_SSE1            = _mm_sub_ps(work->iy_SSE1,jy_SSE);
            dz_SSE1            = _mm_sub_ps(work->iz_SSE1,jz_SSE);
            dx_SSE2            = _mm_sub_ps(work->ix_SSE2,jx_SSE);
            dy_SSE2            = _mm_sub_ps(work->iy_SSE2,jy_SSE);
            dz_SSE2            = _mm_sub_ps(work->iz_SSE2,jz_SSE);
            dx_SSE3            = _mm_sub_ps(work->ix_SSE3,jx_SSE);
            dy_SSE3            = _mm_sub_ps(work->iy_SSE3,jy_SSE);
            dz_SSE3            = _mm_sub_ps(work->iz_SSE3,jz_SSE);
            
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
            
            _mm_store_ps(wco_any_align,wco_any_SSE);
            
            InRange = (wco_any_align[0] != 0 ||
                       wco_any_align[1] != 0 || 
                       wco_any_align[2] != 0 ||
                       wco_any_align[3] != 0);
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
            nbl->cj[nbl->ncj].c    = gridj->cell0 + cj;
            nbl->cj[nbl->ncj].excl = (remove_sub_diag && ci == cj ?
                                      SIMD_INT_MASK_DIAG : SIMD_INT_MASK_ALL);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
#else
    /* No SSE */
    gmx_incons("SSE function called without SSE support");
#endif
}

static void make_subcell_list(const gmx_nbsearch_t nbs,
                              const gmx_nbs_grid_t *gridi,
                              const gmx_nbs_grid_t *gridj,
                              gmx_nblist_t *nbl,
                              int ci,int cj,
                              gmx_bool ci_equals_cj,
                              int stride,const real *x,
                              real rl2,real rbb2)
{
    int  naps;
    int  npair;
    int  sj,si1,si,csj,csj_gl,csi;
    int  sj4_ind,sj_offset;
    unsigned imask;
    gmx_nbl_sj4_t *sj4;
    const real *bb_ci,*x_ci;
    float *d2l,d2;
    int  w;
#define GMX_PRUNE_NBL_CPU_ONE
#ifdef GMX_PRUNE_NBL_CPU_ONE
    int  si_last=-1;
#endif

    d2l = nbl->work->d2;

    bb_ci = nbl->work->bb_ci;
    x_ci  = nbl->work->x_ci;

    naps = gridj->na_tx;

    for(sj=0; sj<gridj->nsubc[cj]; sj++)
    {
        sj4_ind   = (nbl->work->sj_ind >> 2);
        sj_offset = nbl->work->sj_ind - sj4_ind*4;
        sj4       = &nbl->sj4[sj4_ind];
        
        csj = cj*NSUBCELL + sj;

        csj_gl = gridj->cell0*NSUBCELL + csj;

        /* Initialize this j-subcell i-subcell list */
        sj4->sj[sj_offset] = csj_gl;
        imask              = 0;

        if (!nbl->TwoWay && ci_equals_cj)
        {
            si1 = sj + 1;
        }
        else
        {
            si1 = gridi->nsubc[ci];
        }

#ifdef GMX_NBS_BBXXXX
        /* Determine all si1 bb distances in one call with SSE */
        subc_bb_dist2_sse_xxxx(gridj->bb+(csj>>SIMD_WIDTH_2LOG)*NNBSBB_XXXX+(csj & (SIMD_WIDTH-1)),
                               si1,bb_ci,d2l);
#endif

        npair = 0;
        for(si=0; si<si1; si++)
        {
#ifndef GMX_NBS_BBXXXX
            /* Determine the bb distance between csi and csj */
            d2l[si] = subc_bb_dist2(naps,si,bb_ci,csj,gridj->bb);
#endif
            d2 = d2l[si];

/* #define GMX_PRUNE_NBL_CPU */
#ifdef GMX_PRUNE_NBL_CPU
            /* Check if the distance is within the distance where
             * we use only the bounding box distance rbb,
             * or within the cut-off and there is at least one atom pair
             * within the cut-off.
             */
            if (d2 < rbb2 ||
                (d2 < rl2 && nbs->subc_dc(naps,si,x_ci,csj_gl,stride,x,rl2)))
#else
            /* Check if the distance between the two bounding boxes
             * in within the neighborlist cut-off.
             */
            if (d2 < rl2)
#endif
            {
                /* Flag this i-subcell to be taken into account */
                imask |= (1U << (sj_offset*NSUBCELL+si));

#ifdef GMX_PRUNE_NBL_CPU_ONE
                si_last = si;
#endif

                npair++;
            }
        }

#ifdef GMX_PRUNE_NBL_CPU_ONE
        /* If we only found 1 pair, check if any atoms are actually
         * within the cut-off, so we could get rid of it.
         */
        if (npair == 1 && d2l[si_last] >= rbb2)
        {
            if (!nbs->subc_dc(naps,si_last,x_ci,csj_gl,stride,x,rl2))
            {
                imask &= ~(1U << (sj_offset*NSUBCELL+si_last));
                npair--;
            }
        }
#endif

        if (npair > 0)
        {
            /* We have a useful sj entry, close it now */

            /* Set the exclusions for the si== sj entry.
             * Here we don't bother to check if this entry is actually flagged,
             * as it will nearly always be in the list.
             */
            if (ci_equals_cj)
            {
                set_self_and_newton_excls(nbl,sj4_ind,sj_offset,sj);
            }

            /* Copy the sub-cell interaction mask to the list */
            for(w=0; w<NWARP; w++)
            {
                sj4->imei[w].imask |= imask;
            }

            nbl->work->sj_ind++;
            
            /* Keep the count */
            nbl->nsi += npair;

            /* Increase the closing index in i super-cell list */
            nbl->ci[nbl->nci].sj4_ind_end = ((nbl->work->sj_ind+4-1)>>2);
        }
    }
}

static void set_ci_excls(const gmx_nbsearch_t nbs,
                         gmx_nblist_t *nbl,
                         gmx_bool diagRemoved,
                         int na_tx_2log,
                         const gmx_nbl_ci_t *nbl_ci,
                         const t_blocka *excl)
{
    const int *cell;
    int naps;
    int ci;
    int sj_ind_first,sj_ind_last;
    int sj_first,sj_last;
    int ndirect;
    int i,ai,aj,si,eind,ge,se;
    int found,sj_ind_0,sj_ind_1,sj_ind_m;
    int sj_m;
    gmx_bool Found_si;
    int si_ind;
    gmx_nbl_excl_t *nbl_excl;
    int inner_i,inner_e,w;

    cell = nbs->cell;

    naps = nbl->naps;

    if ((nbl->simple && nbl_ci->cj_ind_end == nbl_ci->cj_ind_start) ||
        (!nbl->simple && nbl_ci->sj4_ind_end == nbl_ci->sj4_ind_start))
    {
        /* Empty list */
        return;
    }

    ci = nbl_ci->ci;

    if (nbl->simple)
    {
        sj_ind_first = nbl_ci->cj_ind_start;
        sj_ind_last  = nbl->ncj - 1;

        sj_first = nbl->cj[sj_ind_first].c;
        sj_last  = nbl->cj[sj_ind_last].c;

        /* Determine how many contiguous j-cells we have starting
         * from the first i-cell. This number can be used to directly
         * calculate j-cell indices for excluded atoms.
         */
        ndirect = 0;
        while (sj_ind_first + ndirect <= sj_ind_last &&
               nbl->cj[sj_ind_first+ndirect].c == ci + ndirect)
        {
            ndirect++;
        }
    }
    else
    {
        sj_ind_first = nbl_ci->sj4_ind_start*4;
        sj_ind_last  = nbl->work->sj_ind - 1;

        sj_first = nbl->sj4[nbl_ci->sj4_ind_start].sj[0];
        sj_last  = nbl_sj(nbl,sj_ind_last);

        /* Determine how many contiguous j-sub-cells we have starting
         * from the first i-sub-cell. This number can be used to directly
         * calculate j-sub-cell indices for excluded atoms.
         */
        ndirect = 0;
        while (sj_ind_first + ndirect <= sj_ind_last &&
               nbl_sj(nbl,sj_ind_first+ndirect) == ci*NSUBCELL + ndirect)
        {
            ndirect++;
        }
    }

    /* Loop over the atoms in the i super-cell */
    for(i=0; i<nbl->napc; i++)
    {
        ai = nbs->a[ci*nbl->napc+i];
        if (ai >= 0)
        {
            si  = (i>>na_tx_2log);

            /* Loop over the exclusions for this i-atom */
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
                 * for one-way neighbor lists.
                 */
                if (diagRemoved && ge <= ci*nbl->napc + i)
                {
                    continue;
                }

                se = ge>>na_tx_2log;
                /* Could the sub-cell se be in our list? */
                if (se >= sj_first && se <= sj_last)
                {
                    if (se < sj_first + ndirect)
                    {
                        /* We can calculate sj_ind directly from se */
                        found = sj_ind_first + se - sj_first;
                    }
                    else
                    {
                        /* Search for se using bisection */
                        found = -1;
                        sj_ind_0 = sj_ind_first + ndirect;
                        sj_ind_1 = sj_ind_last + 1;
                        while (found == -1 && sj_ind_0 < sj_ind_1)
                        {
                            sj_ind_m = (sj_ind_0 + sj_ind_1)>>1;

                            if (nbl->simple)
                            {
                                sj_m = nbl->cj[sj_ind_m].c;
                            }
                            else
                            {
                                sj_m = nbl_sj(nbl,sj_ind_m);
                            }

                            if (se == sj_m)
                            {
                                found = sj_ind_m;
                            }
                            else if (se < sj_m)
                            {
                                sj_ind_1 = sj_ind_m;
                            }
                            else
                            {
                                sj_ind_0 = sj_ind_m + 1;
                            }
                        }
                    }

                    if (found >= 0)
                    {
                        inner_i = i  - si*naps;
                        inner_e = ge - se*naps;

                        if (nbl->simple)
                        {
                            nbl->cj[found].excl &= ~(1<<(inner_i*SIMD_WIDTH + inner_e));
                        }
                        else if (nbl_imask0(nbl,found) & (1U << ((found & 3)*NSUBCELL + si)))
                        {
                            w       = (inner_e >> 2);
                            
                            get_nbl_exclusions_1(nbl,found>>2,w,&nbl_excl);
                            
                            nbl_excl->pair[(inner_e & 3)*nbl->naps+inner_i] &=
                                ~(1U << ((found & 3)*NSUBCELL + si));
                        }
                    }
                }
            }
        }
    }
}

static void nb_realloc_ci(gmx_nblist_t *nbl,int n)
{
    nbl->ci_nalloc = over_alloc_small(n);
    nb_realloc_void((void **)&nbl->ci,
                    nbl->nci*sizeof(*nbl->ci),
                    nbl->ci_nalloc*sizeof(*nbl->ci),
                    nbl->alloc,nbl->free);
}

static void new_ci_entry(gmx_nblist_t *nbl,int ci,int shift,int flags,
                         gmx_nbl_work_t *work)
{
    if (nbl->nci + 1 > nbl->ci_nalloc)
    {
        nb_realloc_ci(nbl,nbl->nci+1);
    }
    nbl->ci[nbl->nci].ci            = ci;
    nbl->ci[nbl->nci].shift         = shift;
    if (nbl->simple)
    {
        /* Store the interaction flags along with the shift */
        nbl->ci[nbl->nci].shift    |= flags;
        nbl->ci[nbl->nci].cj_ind_start = nbl->ncj;
        nbl->ci[nbl->nci].cj_ind_end   = nbl->ncj;
    }
    else
    {
        nbl->ci[nbl->nci].sj4_ind_start = nbl->nsj4;
        nbl->ci[nbl->nci].sj4_ind_end   = nbl->nsj4;
    }
}

static void sort_cj_excl(gmx_nbl_cj_t *cj,int ncj,
                         gmx_nbl_work_t *work)
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
        if (cj[j].excl != SIMD_INT_MASK_ALL)
        {
            work->cj[jnew++] = cj[j];
        }
    }
    /* Check if there are exclusions at all or not just the first entry */
    if (!((jnew == 0) || 
          (jnew == 1 && cj[0].excl != SIMD_INT_MASK_ALL)))
    {
        for(j=0; j<ncj; j++)
        {
            if (cj[j].excl == SIMD_INT_MASK_ALL)
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

static void close_ci_entry_simple(gmx_nblist_t *nbl)
{
    int jlen;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    jlen = nbl->ci[nbl->nci].cj_ind_end - nbl->ci[nbl->nci].cj_ind_start;
    if (jlen > 0)
    {
        sort_cj_excl(nbl->cj+nbl->ci[nbl->nci].cj_ind_start,jlen,nbl->work);

        nbl->nci++;
    }
}

static void split_ci_entry(gmx_nblist_t *nbl,
                           int nsp_max_av,gmx_bool progBal,int nc_bal,
                           int thread,int nthread)
{
    int nci_est;
    int nsp_max;
    int sj4_start,sj4_end,j4len,sj4;
    int cci;
    int nsp,nsp_ci,nsp_sj4,nsp_sj4_e,nsp_sj4_p;
    int p;

    /* Estimate the total numbers of ci's of the nblist combined
     * over all threads using the target number of ci's.
     */
    nci_est = nc_bal*thread/nthread + nbl->nci;
    if (progBal)
    {
        /* The first ci blocks should be larger, to avoid overhead.
         * The last ci blocks should be smaller, to improve load balancing.
         */
        nsp_max = max(1,
                      nsp_max_av*nc_bal*3/(2*(nci_est - 1 + nc_bal)));
    }
    else
    {
        nsp_max = nsp_max_av;
    }

    sj4_start = nbl->ci[nbl->nci-1].sj4_ind_start;
    sj4_end   = nbl->ci[nbl->nci-1].sj4_ind_end;
    j4len = sj4_end - sj4_start;

    if (j4len > 1 && j4len*NSUBCELL*4 > nsp_max)
    {
        /* Remove the last ci entry and process the sj4's again */
        nbl->nci -= 1;

        cci       = nbl->nci;
        sj4       = sj4_start;
        nsp       = 0;
        nsp_ci    = 0;
        nsp_sj4_e = 0;
        nsp_sj4   = 0;
        while (sj4 < sj4_end)
        {
            nsp_sj4_p = nsp_sj4;
            nsp_sj4   = 0;
            for(p=0; p<NSUBCELL*4; p++)
            {
                nsp_sj4 += (nbl->sj4[sj4].imei[0].imask >> p) & 1;
            }
            nsp += nsp_sj4;

            if (nsp > nsp_max && nsp > nsp_sj4)
            {
                nbl->ci[cci].sj4_ind_end = sj4;
                cci++;
                nbl->nci++;
                if (nbl->nci+1 > nbl->ci_nalloc)
                {
                    nb_realloc_ci(nbl,nbl->nci+1);
                }
                nbl->ci[cci].ci            = nbl->ci[nbl->nci-1].ci;
                nbl->ci[cci].shift         = nbl->ci[nbl->nci-1].shift;
                nbl->ci[cci].sj4_ind_start = sj4;
                nsp_ci    = nsp - nsp_sj4;
                nsp_sj4_e = nsp_sj4_p;
                nsp       = nsp_sj4;
            }

            sj4++;
        }

        /* Put the remaining sj4's in a new ci entry */
        nbl->ci[cci].sj4_ind_end = sj4_end;

        /* Possibly balance out the last two ci's
         * by moving the last sj4 of the second last ci.
         */
        if (nsp_ci - nsp_sj4_e >= nsp + nsp_sj4_e)
        {
            nbl->ci[cci-1].sj4_ind_end--;
            nbl->ci[cci].sj4_ind_start--;
        }

        cci++;
        nbl->nci++;
    }
}

static void close_ci_entry_supersub(gmx_nblist_t *nbl,
                                    int nsp_max_av,
                                    gmx_bool progBal,int nc_bal,
                                    int thread,int nthread)
{
    int j4len,tlen;
    int nb,b;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    j4len = nbl->ci[nbl->nci].sj4_ind_end - nbl->ci[nbl->nci].sj4_ind_start;
    if (j4len > 0)
    {
        /* We can only have complete blocks of 4 j-entries in a list,
         * so round the count up before closing.
         */
        nbl->nsj4         = ((nbl->work->sj_ind + 4-1) >> 2);
        nbl->work->sj_ind = nbl->nsj4*4;

        nbl->nci++;

        if (nsp_max_av > 0)
        {
            split_ci_entry(nbl,nsp_max_av,progBal,nc_bal,thread,nthread);
        }
    }
}

static void sync_work(gmx_nblist_t *nbl)
{
    nbl->work->sj_ind   = nbl->nsj4*4;
    nbl->work->sj4_init = nbl->nsj4;
}

static void clear_nblist(gmx_nblist_t *nbl)
{
    nbl->nci          = 0;
    nbl->ncj          = 0;
    nbl->nsj4         = 0;
    nbl->nsi          = 0;
    nbl->nexcl        = 1;
}

static void set_icell_bb_simple(const real *bb,int ci,
                                real shx,real shy,real shz,
                                real *bb_ci)
{
    int ia;

    ia = ci*NNBSBB_B;
    bb_ci[0] = bb[ia+0] + shx;
    bb_ci[1] = bb[ia+1] + shy;
    bb_ci[2] = bb[ia+2] + shz;
    bb_ci[4] = bb[ia+4] + shx;
    bb_ci[5] = bb[ia+5] + shy;
    bb_ci[6] = bb[ia+6] + shz;
}

static void set_icell_bb_supersub(const real *bb,int ci,
                                  real shx,real shy,real shz,
                                  real *bb_ci)
{
    int ia,m,i;
    
#ifdef GMX_NBS_BBXXXX
    ia = ci*(NSUBCELL>>SIMD_WIDTH_2LOG)*NNBSBB_XXXX;
    for(m=0; m<(NSUBCELL>>SIMD_WIDTH_2LOG)*NNBSBB_XXXX; m+=NNBSBB_XXXX)
    {
        for(i=0; i<SIMD_WIDTH; i++)
        {
            bb_ci[m+ 0+i] = bb[ia+m+ 0+i] + shx;
            bb_ci[m+ 4+i] = bb[ia+m+ 4+i] + shy;
            bb_ci[m+ 8+i] = bb[ia+m+ 8+i] + shz;
            bb_ci[m+12+i] = bb[ia+m+12+i] + shx;
            bb_ci[m+16+i] = bb[ia+m+16+i] + shy;
            bb_ci[m+20+i] = bb[ia+m+20+i] + shz;
        }
    }
#else
    ia = ci*NSUBCELL*NNBSBB_B;
    for(i=0; i<NSUBCELL*NNBSBB_B; i+=NNBSBB_B)
    {
        bb_ci[i+0] = bb[ia+i+0] + shx;
        bb_ci[i+1] = bb[ia+i+1] + shy;
        bb_ci[i+2] = bb[ia+i+2] + shz;
        bb_ci[i+4] = bb[ia+i+4] + shx;
        bb_ci[i+5] = bb[ia+i+5] + shy;
        bb_ci[i+6] = bb[ia+i+6] + shz;
    }
#endif
}

static void icell_set_x_simple(int ci,
                               real shx,real shy,real shz,
                               int naps,
                               int stride,const real *x,
                               gmx_nbl_work_t *work)
{
    int  ia,i;

    ia = ci*4;

    for(i=0; i<4; i++)
    {
        work->x_ci[i*DIM+0] = x[(ia+i)*stride+XX] + shx;
        work->x_ci[i*DIM+1] = x[(ia+i)*stride+YY] + shy;
        work->x_ci[i*DIM+2] = x[(ia+i)*stride+ZZ] + shz;
    }
}

static void icell_set_x_simple_xxxx(int ci,
                                    real shx,real shy,real shz,
                                    int naps,
                                    int stride,const real *x,
                                    gmx_nbl_work_t *work)
{
    int  ia;

    ia = ci*DIM*SIMD_WIDTH;
#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
    work->ix_SSE0 = _mm_set1_ps(x[ia +  0] + shx);
    work->iy_SSE0 = _mm_set1_ps(x[ia +  4] + shy);
    work->iz_SSE0 = _mm_set1_ps(x[ia +  8] + shz);
    work->ix_SSE1 = _mm_set1_ps(x[ia +  1] + shx);
    work->iy_SSE1 = _mm_set1_ps(x[ia +  5] + shy);
    work->iz_SSE1 = _mm_set1_ps(x[ia +  9] + shz);
    work->ix_SSE2 = _mm_set1_ps(x[ia +  2] + shx);
    work->iy_SSE2 = _mm_set1_ps(x[ia +  6] + shy);
    work->iz_SSE2 = _mm_set1_ps(x[ia + 10] + shz);
    work->ix_SSE3 = _mm_set1_ps(x[ia +  3] + shx);
    work->iy_SSE3 = _mm_set1_ps(x[ia +  7] + shy);
    work->iz_SSE3 = _mm_set1_ps(x[ia + 11] + shz);
#endif
}

static void icell_set_x_supersub(int ci,
                                 real shx,real shy,real shz,
                                 int naps,
                                 int stride,const real *x,
                                 gmx_nbl_work_t *work)
{
    int  ia,i;
    real *x_ci;

    x_ci = work->x_ci;

    ia = ci*NSUBCELL*naps;
    for(i=0; i<NSUBCELL*naps; i++)
    {
        x_ci[i*3 + 0] = x[(ia+i)*stride + 0] + shx;
        x_ci[i*3 + 1] = x[(ia+i)*stride + 1] + shy;
        x_ci[i*3 + 2] = x[(ia+i)*stride + 2] + shz;
    }
}

static void icell_set_x_supersub_sse8(int ci,
                                      real shx,real shy,real shz,
                                      int naps,
                                      int stride,const real *x,
                                      gmx_nbl_work_t *work)
{
    int  si,io,ia,i,j;
    real *x_ci;

    x_ci = work->x_ci;

    for(si=0; si<NSUBCELL; si++)
    {
        for(i=0; i<naps; i+=4)
        {
            io = si*naps + i;
            ia = ci*NSUBCELL*naps + io;
            for(j=0; j<4; j++)
            {
                x_ci[io*3 + j + 0] = x[(ia+j)*stride+0] + shx;
                x_ci[io*3 + j + 4] = x[(ia+j)*stride+1] + shy;
                x_ci[io*3 + j + 8] = x[(ia+j)*stride+2] + shz;
            }
        }
    }
}

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

static int get_nsubpair_max(const gmx_nbsearch_t nbs,
                            int iloc,
                            real rlist,
                            int min_ci_balanced)
{
    const gmx_nbs_grid_t *grid;
    rvec ls;
    real xy_diag2,r_eff_sup,vol_est,nsp_est,nsp_est_nl;
    int  nsubpair_max;

    grid = &nbs->grid[0];

    ls[XX] = (grid->c1[XX] - grid->c0[XX])/(grid->ncx*NSUBCELL_X);
    ls[YY] = (grid->c1[YY] - grid->c0[YY])/(grid->ncy*NSUBCELL_Y);
    ls[ZZ] = (grid->c1[ZZ] - grid->c0[ZZ])*grid->ncx*grid->ncy/(grid->nc*NSUBCELL_Z);

    /* The average squared length of the diagonal of a sub cell */
    xy_diag2 = ls[XX]*ls[XX] + ls[YY]*ls[YY] + ls[ZZ]*ls[ZZ];

    /* The formulas below are a heuristic estimate of the average nsj per si*/
    r_eff_sup = rlist + 0.25*sqrt(xy_diag2);

    if (!nbs->DomDec || nbs->zones->n == 1)
    {
        nsp_est_nl = 0;
    }
    else
    {
        nsp_est_nl =
            sqr(grid->atom_density/grid->na_tx)*
            nonlocal_vol2(nbs->zones,ls,r_eff_sup);
    }

    if (LOCAL_I(iloc))
    {
        /* Sub-cell interacts with itself */
        vol_est  = ls[XX]*ls[YY]*ls[ZZ];
        /* 6/2 rectangular volume on the faces */
        vol_est += xy_diag2*r_eff_sup;
        /* 12/2 quarter pie slices on the edges */
        vol_est += 2*(ls[XX] + ls[YY] + ls[ZZ])*0.25*M_PI*sqr(r_eff_sup);
        /* 4 octants of a sphere */
        vol_est += 0.5*4.0/3.0*M_PI*pow(r_eff_sup,3);

        nsp_est = grid->nsubc_tot*vol_est*grid->atom_density/grid->na_tx;

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

        /* Since the target value is a maximum (this avoid high outlyers,
         * which lead to load imbalance), not average, we get more lists
         * than we ask for (to compensate we need to add NSUBCELL*4/4).
         * But more importantly, the optimal GPU performance moves
         * to lower number of block for very small blocks.
         * To compensate we add the maximum pair count per sj4.
         */
        nsubpair_max += NSUBCELL*4;
    }

    if (debug)
    {
        fprintf(debug,"nbl nsp estimate %.1f, nsubpair_max %d\n",
                nsp_est,nsubpair_max);
    }

    return nsubpair_max;
}

static void print_nblist_ci_cj(FILE *fp,const gmx_nblist_t *nbl)
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
                    nbl->cj[j].c,
                    nbl->cj[j].excl);
        }
    }
}

static void print_nblist_ci_sj(FILE *fp,const gmx_nblist_t *nbl)
{
    int i,j4,j;

    for(i=0; i<nbl->nci; i++)
    {
        fprintf(fp,"ci %4d  shift %2d  nsj4 %2d\n",
                nbl->ci[i].ci,nbl->ci[i].shift,
                nbl->ci[i].sj4_ind_end - nbl->ci[i].sj4_ind_start);

        for(j4=nbl->ci[i].sj4_ind_start; j4<nbl->ci[i].sj4_ind_end; j4++)
        {
            for(j=0; j<4; j++)
            {
                fprintf(fp,"  sj %5d  imask %x\n",
                        nbl->sj4[j4].sj[j],
                        nbl->sj4[j4].imei[0].imask);
            }
        }
    }
}

static void combine_nblists(int nnbl,gmx_nblist_t **nbl,
                            gmx_nblist_t *nblc)
{
    int nci,nsj4,nsi,nexcl;
    int n,i,j4;
    int sj4_offset,si_offset,excl_offset;
    const gmx_nblist_t *nbli;

    nci   = nblc->nci;
    nsj4  = nblc->nsj4;
    nsi   = nblc->nsi;
    nexcl = nblc->nexcl;
    for(i=0; i<nnbl; i++)
    {
        nci   += nbl[i]->nci;
        nsj4  += nbl[i]->nsj4;
        nsi   += nbl[i]->nsi;
        nexcl += nbl[i]->nexcl;
    }

    if (nci > nblc->ci_nalloc)
    {
        nb_realloc_ci(nblc,nci);
    }
    if (nsj4 > nblc->sj4_nalloc)
    {
        nblc->sj4_nalloc = over_alloc_small(nsj4);
        nb_realloc_void((void **)&nblc->sj4,
                        nblc->nsj4*sizeof(*nblc->sj4),
                        nblc->sj4_nalloc*sizeof(*nblc->sj4),
                        nblc->alloc,nblc->free);
    }
    if (nexcl > nblc->excl_nalloc)
    {
        nblc->excl_nalloc = over_alloc_small(nexcl);
        nb_realloc_void((void **)&nblc->excl,
                        nblc->nexcl*sizeof(*nblc->excl),
                        nblc->excl_nalloc*sizeof(*nblc->excl),
                        nblc->alloc,nblc->free);
    }

    for(n=0; n<nnbl; n++)
    {
        sj4_offset  = nblc->nsj4;
        si_offset   = nblc->nsi;
        excl_offset = nblc->nexcl;

        nbli = nbl[n];

        /* We could instead omp prallelizing the two loops below
         * parallelize the copy of integral parts of the nblist.
         * However this requires a lot more bookkeeping and does not
         * lead to a performance improvement.
         */
        /* The ci list copy is probably not worth parallelizing */
        for(i=0; i<nbli->nci; i++)
        {
            nblc->ci[nblc->nci]                = nbli->ci[i];
            nblc->ci[nblc->nci].sj4_ind_start += sj4_offset;
            nblc->ci[nblc->nci].sj4_ind_end   += sj4_offset;
            nblc->nci++;
        }

#pragma omp parallel
        {
#pragma omp for schedule(static) nowait
            for(j4=0; j4<nbli->nsj4; j4++)
            {
                nblc->sj4[nblc->nsj4+j4] = nbli->sj4[j4];
                nblc->sj4[nblc->nsj4+j4].imei[0].excl_ind += excl_offset;
                nblc->sj4[nblc->nsj4+j4].imei[1].excl_ind += excl_offset;
            }

#pragma omp for schedule(static)
            for(j4=0; j4<nbli->nexcl; j4++)
            {
                nblc->excl[nblc->nexcl+j4] = nbli->excl[j4];
            }
        }

        nblc->nsj4  += nbli->nsj4;
        nblc->nsi   += nbli->nsi;
        nblc->nexcl += nbli->nexcl;
    }
}

static gmx_bool next_ci(const gmx_nbs_grid_t *grid,
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

static real boundingbox_only_distance2(const gmx_nbs_grid_t *gridi,
                                       const gmx_nbs_grid_t *gridj,
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

    bbx = 0.5*(gridi->sx + gridj->sx);
    bby = 0.5*(gridi->sy + gridj->sy);
    if (!simple)
    {
        bbx /= NSUBCELL_X;
        bby /= NSUBCELL_Y;
    }

    return sqr(max(0,rlist - 0.5*sqrt(bbx*bbx + bby*bby)));
}

static void gmx_nbsearch_make_nblist_part(const gmx_nbsearch_t nbs,
                                          const gmx_nbs_grid_t *gridi,
                                          const gmx_nbs_grid_t *gridj,
                                          gmx_nbs_work_t *work,
                                          const gmx_nb_atomdata_t *nbat,
                                          const t_blocka *excl,
                                          real rlist,
                                          int nsubpair_max,
                                          gmx_bool progBal,
                                          int min_ci_balanced,
                                          int th,int nth,
                                          gmx_nblist_t *nbl)
{
    matrix box;
    real rl2,rbb2;
    int  d;
    int  ci_block,ci_b,ci,ci_x,ci_y,ci_xy,cj;
    ivec shp;
    int  tx,ty,tz;
    int  shift;
    gmx_bool bMakeList;
    real shx,shy,shz;
    int  conv_i,cell0_i;
    const gmx_bool *use_c;
    const real *bb_i,*bbcz_i,*bbcz_j;
    const int *flags_i;
    real bx0,bx1,by0,by1,bz0,bz1;
    real bz1_frac;
    real d2cx,d2z,d2z_cx,d2z_cy,d2zx,d2zxy,d2xy;
    int  cxf,cxl,cyf,cyf_x,cyl;
    int  cx,cy;
    int  c0,c1,cs,cf,cl;
    int  ncpc;

    nbs_cycle_start(&work->cc[enbsCCsearch]);

    sync_work(nbl);

    /* Maybe we should not set the type of neighborlist at search time,
     * but rather in gmx_nblist_init?
     */
    nbl->simple = gridj->simple;

    nbl->napc = gridj->na_stx;
    nbl->naps = gridj->na_tx;

    /* Currently this code only makes two-way lists */
    nbl->TwoWay = FALSE;

    nbl->rlist  = rlist;

    copy_mat(nbs->box,box);

    rl2 = nbl->rlist*nbl->rlist;

    rbb2 = boundingbox_only_distance2(gridi,gridj,nbl->rlist,nbl->simple);

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

    if (nbl->simple && !gridi->simple)
    {
        conv_i  = gridi->na_stx/gridj->na_stx;
        use_c   = gridi->use_simple;
        bb_i    = gridi->bb_simple;
        bbcz_i  = gridi->bbcz_simple;
        flags_i = gridi->flags_simple;
    }
    else
    {
        conv_i  = 1;
        use_c   = NULL;
        bb_i    = gridi->bb;
        bbcz_i  = gridi->bbcz;
        flags_i = gridi->flags;
    }
    cell0_i = gridi->cell0*conv_i;

    bbcz_j = gridj->bbcz;

    if (conv_i == 1)
    {
        /* Set the block size as 5/11/ntask times the average number of cells
         * in a y,z slab. This should ensure a quite uniform distribution
         * of the grid parts of the different thread along all three grid
         * zone boundaries with 3D domain decomposition. At the same time
         * the blocks will not become too small.
         */
        ci_block = (gridi->nc*5)/(11*gridi->ncx*nth);

        /* Ensure the blocks are not too small: avoids cache invalidation */
        if (ci_block*gridi->na_stx < 16)
        {
            ci_block = (16 + gridi->na_stx - 1)/gridi->na_stx;
        }

        /* Without domain decomposition
         * or with less than 3 blocks per task, divide in nth blocks.
         */
        if (!nbs->DomDec || ci_block*3*nth > gridi->nc)
        {
            ci_block = (gridi->nc + nth - 1)/nth;
        }
    }
    else
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

    ncpc = 0;

    ci_b = -1;
    ci   = th*ci_block - 1;
    ci_x = 0;
    ci_y = 0;
    while (next_ci(gridi,conv_i,nth,ci_block,&ci_x,&ci_y,&ci_b,&ci))
    {
        if (conv_i > 1 && !use_c[ci])
        {
            continue;
        }

        d2cx = 0;
        if (gridj != gridi && shp[XX] == 0)
        {
            if (nbl->simple)
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

                if (nbl->simple)
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

#ifdef NSBOX_SHIFT_BACKWARD
                    if (gridi == gridj && shift > CENTRAL)
                    {
                        continue;
                    }
#endif

                    shx = tx*box[XX][XX] + ty*box[YY][XX] + tz*box[ZZ][XX];

                    if (nbl->simple)
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

                    new_ci_entry(nbl,cell0_i+ci,shift,flags_i[ci],
                                 nbl->work);

#ifndef NSBOX_SHIFT_BACKWARD
                    if (!nbl->TwoWay && cxf < ci_x)
#else
                    if (!nbl->TwoWay && shift == CENTRAL && gridi == gridj &&
                        cxf < ci_x)
#endif
                    {
                        /* Leave the pairs with i > j.
                         * x is the major index, so skip half of it.
                         */
                        cxf = ci_x;
                    }

                    if (nbl->simple)
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
                                     gridi->na_tx,nbat->xstride,nbat->x,nbl->work);

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

#ifndef NSBOX_SHIFT_BACKWARD
                        if (!nbl->TwoWay && gridi == gridj &&
                            cx == 0 && cyf < ci_y)
#else
                        if (!nbl->TwoWay && gridi == gridj &&
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
#ifdef NSBOX_SHIFT_BACKWARD
                            if (!nbl->TwoWay && gridi == gridj &&
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

#ifdef NSBOX_REFCODE
                                {
                                    /* Simple reference code */
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
#ifndef NSBOX_SHIFT_BACKWARD
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
                                    if (nbl->simple)
                                    {
                                        check_subcell_list_space_simple(nbl,cl-cf+1);

                                        if (nbat->XFormat == nbatXXXX)
                                        {
                                            make_subcell_list_simple_xxxx(gridj,
                                                                          nbl,ci,cf,cl,
                                                                          (gridi == gridj && shift == CENTRAL),
                                                                          nbat->x,
                                                                          rl2,rbb2);
                                        }
                                        else
                                        {
                                            make_subcell_list_simple(gridj,
                                                                     nbl,ci,cf,cl,
                                                                     (gridi == gridj && shift == CENTRAL),
                                                                     nbat->x,
                                                                     rl2,rbb2);
                                        }
                                    }
                                    else
                                    {
                                        check_subcell_list_space_supersub(nbl,cl-cf+1);

                                        for(cj=cf; cj<=cl; cj++)
                                        {
                                            make_subcell_list(nbs,gridi,gridj,
                                                              nbl,ci,cj,
                                                              (gridi == gridj && shift == CENTRAL && ci == cj),
                                                              nbat->xstride,nbat->x,
                                                              rl2,rbb2);
                                        }
                                    }
                                    ncpc += cl - cf + 1;
                                }
                            }
                        }  
                    }

                    /* Set the exclusions for this ci list */
                    set_ci_excls(nbs,
                                 nbl,
                                 !nbl->TwoWay && shift == CENTRAL && gridi == gridj,
                                 gridj->na_tx_2log,
                                 &(nbl->ci[nbl->nci]),excl);

                    /* Close this ci list */
                    if (nbl->simple)
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
    }

    nbs_cycle_stop(&work->cc[enbsCCsearch]);

    if (debug)
    {
        fprintf(debug,"ncpc %s %d\n",gridi==gridj ? "local" : "non-local",ncpc);

        if (nbl->simple)
        {
            print_nblist_statistics_simple(debug,nbl,nbs,rlist);
        }
        else
        {
            print_nblist_statistics_supersub(debug,nbl,nbs,rlist);
        }

    }
}

void gmx_nbsearch_make_nblist(const gmx_nbsearch_t nbs,
                              const gmx_nb_atomdata_t *nbat,
                              const t_blocka *excl,
                              real rlist,
                              int min_ci_balanced,
                              gmx_nbl_lists_t *nbl_list,
                              int iloc)
{
    const gmx_nbs_grid_t *gridi,*gridj;
    int nzi,zi,zj0,zj1,zj;
    int nsubpair_max;
    int nth,th;
    int nnbl;
    gmx_nblist_t **nbl;
    gmx_bool CombineNBLists;

    nnbl            = nbl_list->nnbl;
    nbl             = nbl_list->nbl;
    CombineNBLists  = nbl_list->combined;

    if (debug)
    {
        fprintf(debug,"ns making %d nblists\n", nnbl);
    }

    if (nbl_list->simple)
    {
        if (nbat->XFormat == nbatXXXX)
        {
            nbs->icell_set_x = icell_set_x_simple_xxxx;
        }
        else
        {
            nbs->icell_set_x = icell_set_x_simple;
        }
    }
    else
    {
#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
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

    if (!nbl_list->simple && min_ci_balanced > 0)
    {
        nsubpair_max = get_nsubpair_max(nbs,iloc,rlist,min_ci_balanced);
    }
    else
    {
        nsubpair_max = 0;
    }

    /* Clear all the neighbor lists */
    for(th=0; th<nnbl; th++)
    {
        clear_nblist(nbl[th]);
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

#pragma omp parallel for schedule(static)
            for(th=0; th<nnbl; th++)
            {
                if (CombineNBLists && th > 0)
                {
                    clear_nblist(nbl[th]);
                }

                /* Divide the i super cell equally over the nblists */
                gmx_nbsearch_make_nblist_part(nbs,gridi,gridj,
                                              &nbs->work[th],nbat,excl,
                                              rlist,
                                              nsubpair_max,
                                              (LOCAL_I(iloc) || nbs->zones->n <= 2),
                                              min_ci_balanced,
                                              th,nnbl,
                                              nbl[th]);
            }
            nbs_cycle_stop(&nbs->cc[enbsCCsearch]);

            if (CombineNBLists && nnbl > 1)
            {
                nbs_cycle_start(&nbs->cc[enbsCCcombine]);

                combine_nblists(nnbl-1,nbl+1,nbl[0]);

                nbs_cycle_stop(&nbs->cc[enbsCCcombine]);
            }
        }
    }

    /*
    print_supersub_nsp("nsubpair",nbl[0],iloc);
    */

    if (nbs->print_cycles &&
        nbs->cc[enbsCCgrid].count > 0 &&
        nbs->cc[enbsCCgrid].count % 100 == 0)
    {
        nbs_cycle_print(stderr,nbs);
    }

    if (debug && (CombineNBLists && nnbl > 1))
    {
        if (nbl[0]->simple)
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
        if (nbl[0]->simple)
        {
            print_nblist_ci_cj(debug,nbl[0]);
        }
        else
        {
            print_nblist_ci_sj(debug,nbl[0]);
        }   
    }
}

static void gmx_nb_atomdata_output_init(gmx_nb_atomdata_output_t *out,
                                        gmx_bool simple,
                                        int nenergrp,
                                        gmx_nbat_alloc_t *ma)
{
    out->f = NULL;
    ma((void **)&out->fshift,SHIFTS*DIM*sizeof(*out->fshift));
    out->nV = nenergrp*nenergrp;
    ma((void **)&out->Vvdw,out->nV*sizeof(*out->Vvdw));
    ma((void **)&out->Vc  ,out->nV*sizeof(*out->Vc  ));
    if (simple)
    {
        out->nVS = nenergrp*nenergrp*nenergrp*nenergrp*nenergrp*4;
        ma((void **)&out->VSvdw,out->nVS*sizeof(*out->VSvdw));
        ma((void **)&out->VSc  ,out->nVS*sizeof(*out->VSc  ));
    }
    else
    {
        out->nVS = 0;
    }
}

static void set_combination_rule_data(gmx_nb_atomdata_t *nbat)
{
    int  nt,i,j;
    real c6,c12;

    nt = nbat->ntype;

    switch (nbat->comb_rule)
    {
    case  ljcrGEOM:
        nbat->comb_rule = ljcrGEOM;
        
        for(i=0; i<nt; i++)
        {
            /* Copy the diagonal from the nbfp matrix */
            nbat->nbfp_comb[i*2  ] = sqrt(nbat->nbfp[(i*nt+i)*2  ]);
            nbat->nbfp_comb[i*2+1] = sqrt(nbat->nbfp[(i*nt+i)*2+1]);
        }
        break;
    case ljcrLB:
        for(i=0; i<nt; i++)
        {
            /* Get 6*C6 and 12*C12 from the diagonal of the nbfp matrix */
            c6  = nbat->nbfp[(i*nt+i)*2  ];
            c12 = nbat->nbfp[(i*nt+i)*2+1];
            if (c6 > 0 && c12 > 0)
            {
                /* We store 0.5*2^1/6*sigma and sqrt(4*3*eps),
                 * so we get 6*C6 and 12*C12 after combining.
                 */
                nbat->nbfp_comb[i*2  ] = 0.5*pow(c12/c6,1.0/6.0);
                nbat->nbfp_comb[i*2+1] = sqrt(c6*c6/c12);
            }
            else
            {
                nbat->nbfp_comb[i*2  ] = 0;
                nbat->nbfp_comb[i*2+1] = 0;
            }
        }
        break;
    case ljcrNONE:
        nbat->alloc((void **)&nbat->nbfp_s4,nt*nt*4*sizeof(*nbat->nbfp_s4));
        for(i=0; i<nt; i++)
        {
            for(j=0; j<nt; j++)
            {
                nbat->nbfp_s4[(i*nt+j)*4+0] = nbat->nbfp[(i*nt+j)*2+0];
                nbat->nbfp_s4[(i*nt+j)*4+1] = nbat->nbfp[(i*nt+j)*2+1];
                nbat->nbfp_s4[(i*nt+j)*4+2] = 0;
                nbat->nbfp_s4[(i*nt+j)*4+3] = 0;
            }
        }
        break;
    default:
        gmx_incons("Unknown combination rule");
        break;
    }
}

void gmx_nb_atomdata_init(FILE *fp,
                          gmx_nb_atomdata_t *nbat,
                          gmx_bool simple,
                          gmx_bool XFormatXXXX,
                          int ntype,const real *nbfp,
                          int n_energygroups,
                          int nout,
                          gmx_nbat_alloc_t *alloc,
                          gmx_nbat_free_t  *free)
{
    int  i,j;
    real c6,c12,tol;
    char *ptr;
    gmx_bool bCombGeom,bCombLB;

    if (alloc == NULL)
    {
        nbat->alloc = nblist_alloc_aligned;
    }
    else
    {
        nbat->alloc = alloc;
    }
    if (free == NULL)
    {
        nbat->free = nblist_free_aligned;
    }
    else
    {
        nbat->free = free;
    }

    if (debug)
    {
        fprintf(debug,"There are %d atom types in the system, adding one for gmx_nb_atomdata_t\n",ntype);
    }
    nbat->ntype = ntype + 1;
    nbat->alloc((void **)&nbat->nbfp,
                nbat->ntype*nbat->ntype*2*sizeof(*nbat->nbfp));
    nbat->alloc((void **)&nbat->nbfp_comb,nbat->ntype*2*sizeof(*nbat->nbfp_comb));

    /* A tolerance of 1e-5 seems reasonable for (possibly hand-typed)
     * force-field floating point parameters.
     */
    tol = 1e-5;
    ptr = getenv("GMX_LJCOMB_TOL");
    if (ptr != NULL)
    {
        double dbl;

        sscanf(ptr,"%lf",&dbl);
        tol = dbl;
    }
    bCombGeom = TRUE;
    bCombLB   = TRUE;

    /* Temporarily fill nbat->nbfp_comb with sigma and epsilon
     * to check for the LB rule.
     */
    for(i=0; i<ntype; i++)
    {
        c6  = nbfp[(i*ntype+i)*2  ];
        c12 = nbfp[(i*ntype+i)*2+1];
        if (c6 > 0 && c12 > 0)
        {
            nbat->nbfp_comb[i*2  ] = pow(c12/c6,1.0/6.0);
            nbat->nbfp_comb[i*2+1] = 0.25*c6*c6/c12;
        }
        else if (c6 == 0 && c12 == 0)
        {
            nbat->nbfp_comb[i*2  ] = 0;
            nbat->nbfp_comb[i*2+1] = 0;
        }
        else
        {
            /* Can not use LB rule with only dispersion or repulsion */
            bCombLB = FALSE;
        }
    }

    for(i=0; i<nbat->ntype; i++)
    {
        for(j=0; j<nbat->ntype; j++)
        {
            if (i < ntype && j < ntype)
            {
                /* We store the prefactor in the derivative of the potential
                 * in the parameter to avoid multiplications in the inner loop.
                 */
                c6  = nbfp[(i*ntype+j)*2  ];
                c12 = nbfp[(i*ntype+j)*2+1];
                nbat->nbfp[(i*nbat->ntype+j)*2  ] =  6.0*c6;
                nbat->nbfp[(i*nbat->ntype+j)*2+1] = 12.0*c12;

                bCombGeom = bCombGeom &&
                    gmx_within_tol(c6*c6  ,nbfp[(i*ntype+i)*2  ]*nbfp[(j*ntype+j)*2  ],tol) &&
                    gmx_within_tol(c12*c12,nbfp[(i*ntype+i)*2+1]*nbfp[(j*ntype+j)*2+1],tol);

                bCombLB = bCombLB &&
                    ((c6 == 0 && c12 == 0 &&
                      (nbat->nbfp_comb[i*2+1] == 0 || nbat->nbfp_comb[j*2+1] == 0)) ||
                     (c6 > 0 && c12 > 0 &&
                      gmx_within_tol(pow(c12/c6,1.0/6.0),0.5*(nbat->nbfp_comb[i*2]+nbat->nbfp_comb[j*2]),tol) &&
                      gmx_within_tol(0.25*c6*c6/c12,sqrt(nbat->nbfp_comb[i*2+1]*nbat->nbfp_comb[j*2+1]),tol)));
            }
            else
            {
                /* Add zero parameters for the additional dummy atom type */
                nbat->nbfp[(i*nbat->ntype+j)*2  ] = 0;
                nbat->nbfp[(i*nbat->ntype+j)*2+1] = 0;
            }
        }
    }
    if (debug)
    {
        fprintf(debug,"Combination rules: geometric %d Lorentz-Berthelot %d\n",
                bCombGeom,bCombLB);
    }

    /* This should be switched by the kernel type, not nblist type */
    if (simple)
    {
        /* We prefer the geometic combination rule,
         * as that give a slightly faster kernel than the LB rule.
         */
        if (bCombGeom)
        {
            nbat->comb_rule = ljcrGEOM;
        }
        else if (bCombLB)
        {
            nbat->comb_rule = ljcrLB;
        }
        else
        {
            nbat->comb_rule = ljcrNONE;

            nbat->free(nbat->nbfp_comb);           
        }

        if (fp)
        {
            if (nbat->comb_rule == ljcrNONE)
            {
                fprintf(fp,"Using full Lennard-Jones parameter combination matrix\n\n");
            }
            else
            {
                fprintf(fp,"Using %s Lennard-Jones combination rule\n\n",
                        nbat->comb_rule==ljcrGEOM ? "geometric" : "Lorentz-Berthelot");
            }
        }

        set_combination_rule_data(nbat);
    }
    else
    {
        nbat->comb_rule = ljcrNONE;

        nbat->free(nbat->nbfp_comb);
    }

    nbat->natoms  = 0;
    nbat->type    = NULL;
    nbat->lj_comb = NULL;
    if (simple)
    {
        nbat->XFormat = (XFormatXXXX ? nbatXXXX : nbatXYZ);
    }
    else
    {
        nbat->XFormat = nbatXYZQ;
    }
    nbat->q       = NULL;
    nbat->nenergrp = n_energygroups;
    if (!simple)
    {
        /* Energy groups not supported yet for super-sub lists */
        nbat->nenergrp = 1;
    }
    nbat->energrp = NULL;
    nbat->alloc((void **)&nbat->shift_vec,SHIFTS*sizeof(*nbat->shift_vec));
    nbat->xstride = (nbat->XFormat == nbatXYZQ ? 4 : 3);
    nbat->x       = NULL;
    nbat->nout    = nout;
    snew(nbat->out,nbat->nout);
    nbat->nalloc  = 0;
    for(i=0; i<nbat->nout; i++)
    {
        gmx_nb_atomdata_output_init(&nbat->out[i],
                                    simple,nbat->nenergrp,nbat->alloc);
    }
}

static void copy_lj_to_nbat_lj_comb(const real *ljparam_type,
                                    const int *type,int na,
                                    real *ljparam_at)
{
    int is,k,i;

    /* The LJ params follow the combination rule:
     * copy the params for the type array to the atom array.
     */
    for(is=0; is<na; is+=SIMD_WIDTH)
    {
        for(k=0; k<SIMD_WIDTH; k++)
        {
            i = is + k;
            ljparam_at[is*2           +k] = ljparam_type[type[i]*2  ];
            ljparam_at[is*2+SIMD_WIDTH+k] = ljparam_type[type[i]*2+1];
        }
    }
}

static void gmx_nb_atomdata_set_atomtypes(gmx_nb_atomdata_t *nbat,
                                          int ngrid,
                                          const gmx_nbsearch_t nbs,
                                          const int *type)
{
    int g,i,ncz,ash;
    const gmx_nbs_grid_t *grid;

    for(g=0; g<ngrid; g++)
    {
        grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for(i=0; i<grid->ncx*grid->ncy; i++)
        {
            ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
            ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_stx;

            copy_int_to_nbat_int(nbs->a+ash,grid->cxy_na[i],ncz*grid->na_stx,
                                 type,nbat->ntype-1,nbat->type+ash);

            if (nbat->comb_rule != ljcrNONE)
            {
                copy_lj_to_nbat_lj_comb(nbat->nbfp_comb,
                                        nbat->type+ash,ncz*grid->na_stx,
                                        nbat->lj_comb+ash*2);
            }
        }
    }
}

static void gmx_nb_atomdata_set_charges(gmx_nb_atomdata_t *nbat,
                                        int ngrid,
                                        const gmx_nbsearch_t nbs,
                                        const real *charge)
{
    int  g,cxy,ncz,ash,na,na_round,i,j;
    real *q;
    const gmx_nbs_grid_t *grid;

    for(g=0; g<ngrid; g++)
    {
        grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for(cxy=0; cxy<grid->ncx*grid->ncy; cxy++)
        {
            ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_stx;
            na  = grid->cxy_na[cxy];
            na_round = (grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy])*grid->na_stx;

            if (nbat->XFormat == nbatXYZQ)
            {
                q = nbat->x + ash*nbat->xstride + 3;
                for(i=0; i<na; i++)
                {
                    *q = charge[nbs->a[ash+i]];
                    q += 4;
                }
                /* Complete the partially filled last cell with zeros */
                for(; i<na_round; i++)
                {
                    *q = 0;
                    q += 4;
                }
            }
            else
            {
                q = nbat->q + ash;
                for(i=0; i<na; i++)
                {
                    *q = charge[nbs->a[ash+i]];
                    q++;
                }
                /* Complete the partially filled last cell with zeros */
                for(; i<na_round; i++)
                {
                    *q = 0;
                    q++;
                }
            }
        }
    }
}

static void copy_egp_to_nbat_egps(const int *a,int na,int na_round,
                                  int naps,int comb_fac,
                                  const int *in,int *innb)
{
    int i,j,sa,at;
    int comb,fac;

    j = 0;
    for(i=0; i<na; i+=naps)
    {
        /* Store naps energy groups number into one int */
        comb = 0;
        fac  = 1;
        for(sa=0; sa<naps; sa++)
        {
            at = a[i+sa];
            if (at >= 0)
            {
                comb += GET_CGINFO_GID(in[at])*fac;
            }
            fac *= comb_fac;
        }
        innb[j++] = comb;
    }
    /* Complete the partially filled last cell with fill */
    for(; i<na_round; i+=naps)
    {
        innb[j++] = 0;
    }
}

static void gmx_nb_atomdata_set_energygroups(gmx_nb_atomdata_t *nbat,
                                             int ngrid,
                                             const gmx_nbsearch_t nbs,
                                             const int *atinfo)
{
    int g,i,ncz,ash;
    const gmx_nbs_grid_t *grid;
    int comb_fac;

    if (nbat->XFormat == nbatXXXX)
    {
        /* Minimize the maximum combined energy group number */
        comb_fac = nbat->nenergrp;
    }
    else
    {
        /* Easy access through bit operations by using a power of 2 */
        comb_fac = (1<<8);
    }

    for(g=0; g<ngrid; g++)
    {
        grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for(i=0; i<grid->ncx*grid->ncy; i++)
        {
            ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
            ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_stx;

            copy_egp_to_nbat_egps(nbs->a+ash,grid->cxy_na[i],ncz*grid->na_stx,
                                  nbat->naps,comb_fac,
                                  atinfo,nbat->energrp+(ash>>grid->na_tx_2log));
        }
    }
}

void gmx_nb_atomdata_set(gmx_nb_atomdata_t *nbat,
                         int aloc,
                         const gmx_nbsearch_t nbs,
                         const t_mdatoms *mdatoms,
                         const int *atinfo)
{
    int ngrid;

    if (aloc == eatLocal)
    {
        ngrid = 1;
    }
    else
    {
        ngrid = nbs->ngrid;
    }

    gmx_nb_atomdata_set_atomtypes(nbat,ngrid,nbs,mdatoms->typeA);

    gmx_nb_atomdata_set_charges(nbat,ngrid,nbs,mdatoms->chargeA);

    if (nbat->nenergrp > 1)
    {
        gmx_nb_atomdata_set_energygroups(nbat,ngrid,nbs,atinfo);
    }
}

void gmx_nb_atomdata_copy_shiftvec(gmx_bool dynamic_box,
                                   rvec *shift_vec,
                                   gmx_nb_atomdata_t *nbat)
{
    int i;

    nbat->dynamic_box = dynamic_box;
    for(i=0; i<SHIFTS; i++)
    {
        copy_rvec(shift_vec[i],nbat->shift_vec[i]);
    }
}

void gmx_nb_atomdata_copy_x_to_nbat_x(const gmx_nbsearch_t nbs,
                                      int aloc,
                                      gmx_bool FillLocal,
                                      rvec *x,
                                      gmx_nb_atomdata_t *nbat)
{
    int g0=0,g1=0,g,cxy,na,ash,na_fill;
    const gmx_nbs_grid_t *grid;

    switch (aloc)
    {
    case eatAll:
        g0 = 0;
        g1 = nbs->ngrid;
        break;
    case eatLocal:
        g0 = 0;
        g1 = 1;
        break;
    case eatNonlocal:
        g0 = 1;
        g1 = nbs->ngrid;
        break;
    }

    if (FillLocal)
    {
        nbat->natoms_local = nbs->grid[0].nc*nbs->grid[0].na_stx;
    }

    for(g=g0; g<g1; g++)
    {
        grid = &nbs->grid[g];

#pragma omp parallel for schedule(static) private(na,ash)
        for(cxy=0; cxy<grid->ncx*grid->ncy; cxy++)
        {
            na  = grid->cxy_na[cxy];
            ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_stx;

            if (g == 0 && FillLocal)
            {
                na_fill =
                    (grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy])*grid->na_stx;
            }
            else
            {
                /* We fill only the real particle locations.
                 * We assume the filling entries at the end have been
                 * properly set before during ns.
                 */
                na_fill = na;
            }
            copy_rvec_to_nbat_real(nbs->a+ash,na,na_fill,x,
                                   nbat->XFormat,nbat->x+ash*nbat->xstride,
                                   0,0,0);
        }
    }
}

static void gmx_nb_atomdata_add_nbat_f_to_f_part(const gmx_nbsearch_t nbs,
                                                 const gmx_nb_atomdata_t *nbat,
                                                 gmx_nb_atomdata_output_t *out,
                                                 int nfa,
                                                 int a0,int a1,
                                                 rvec *f)
{
    int  a,i,fa;
    const int  *cell;
    const real *fnb;

    cell = nbs->cell;

    /* Loop over all columns and copy and fill */
    switch (nbat->XFormat)
    {
    case nbatXYZ:
    case nbatXYZQ:
        if (nfa == 1)
        {
            fnb = out[0].f;

            for(a=a0; a<a1; a++)
            {
                i = cell[a]*nbat->xstride;
                
                f[a][XX] += fnb[i];
                f[a][YY] += fnb[i+1];
                f[a][ZZ] += fnb[i+2];
            }
        }
        else
        {
            for(a=a0; a<a1; a++)
            {
                i = cell[a]*nbat->xstride;
                
                for(fa=0; fa<nfa; fa++)
                {
                    f[a][XX] += out[fa].f[i];
                    f[a][YY] += out[fa].f[i+1];
                    f[a][ZZ] += out[fa].f[i+2];
                } 
            }
        }
        break;
    case nbatXXXX:
        if (nfa == 1)
        {
            fnb = out[0].f;

            for(a=a0; a<a1; a++)
            {
                i = 12*(cell[a] >> SIMD_WIDTH_2LOG) +
                    (cell[a] & (SIMD_WIDTH-1));

                f[a][XX] += fnb[i];
                f[a][YY] += fnb[i+4];
                f[a][ZZ] += fnb[i+8];
            }
        }
        else
        {
            for(a=a0; a<a1; a++)
            {
                i = 12*(cell[a] >> SIMD_WIDTH_2LOG) +
                    (cell[a] & (SIMD_WIDTH-1));
                
                for(fa=0; fa<nfa; fa++)
                {
                    f[a][XX] += out[fa].f[i];
                    f[a][YY] += out[fa].f[i+4];
                    f[a][ZZ] += out[fa].f[i+8];
                }
            }
        }
        break;
    }
}

void gmx_nb_atomdata_add_nbat_f_to_f(const gmx_nbsearch_t nbs,
                                     int aloc,
                                     const gmx_nb_atomdata_t *nbat,
                                     rvec *f)
{
    int a0,na;
    int nth,th;

    nbs_cycle_start(&nbs->cc[enbsCCreducef]);

    switch (aloc)
    {
    case eatAll:
        a0 = 0;
        na = nbs->natoms_nonlocal;
        break;
    case eatLocal:
        a0 = 0;
        na = nbs->natoms_local;
        break;
    case eatNonlocal:
        a0 = nbs->natoms_local;
        na = nbs->natoms_nonlocal - nbs->natoms_local;
        break;
    }

    nth = omp_get_max_threads();
#pragma omp parallel for schedule(static)
    for(th=0; th<nth; th++)
    {
        gmx_nb_atomdata_add_nbat_f_to_f_part(nbs,nbat,
                                             nbat->out,
                                             nbat->nout,
                                             a0+((th+0)*na)/nth,
                                             a0+((th+1)*na)/nth,
                                             f);
    }

    nbs_cycle_stop(&nbs->cc[enbsCCreducef]);
}

void gmx_nb_atomdata_add_nbat_fshift_to_fshift(const gmx_nb_atomdata_t *nbat,
                                               rvec *fshift)
{
    const gmx_nb_atomdata_output_t *out;
    int  nth,th;
    int  s;
    rvec sum;

    out = nbat->out;

    nth = omp_get_max_threads();
    
    for(s=0; s<SHIFTS; s++)
    {
        clear_rvec(sum);
        for(th=0; th<nth; th++)
        {
            sum[XX] += out[th].fshift[s*DIM+XX];
            sum[YY] += out[th].fshift[s*DIM+YY];
            sum[ZZ] += out[th].fshift[s*DIM+ZZ];
        }
        rvec_inc(fshift[s],sum);
    }
}
