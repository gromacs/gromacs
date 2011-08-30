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

    gmx_cache_protect_t cp1;
} gmx_nbl_work_t;

typedef void
gmx_supcell_set_i_x_t(const gmx_nbsearch_t nbs,int ci,
                      real shx,real shy,real shz,
                      int stride,const real *x,
                      gmx_nbl_work_t *work);

static gmx_supcell_set_i_x_t supcell_set_i_x;
static gmx_supcell_set_i_x_t supcell_set_i_x_sse8;

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

enum { enbsCCgrid, enbsCCsearch, enbsCCcombine, enbsCCnr };

typedef struct {
    gmx_cache_protect_t cp0;

    int *cxy_na;

    int  *sort_work;
    int  sort_work_nalloc;

    float *bb_tmp;

    gmx_nbs_cc_t cc[enbsCCnr];

    gmx_cache_protect_t cp1;
} gmx_nbs_work_t;

typedef struct gmx_nbsearch {
    int  ePBC;
    matrix box;
    real atom_density;

    int  naps;  /* Number of atoms in the inner loop / sub cell */
    int  napc;  /* Number of atoms in the super cell            */
    int  naps2log;

    int  ncx;
    int  ncy;
    int  nc;

    real sx;
    real sy;
    real inv_sx;
    real inv_sy;

    int  *cxy_na;
    int  *cxy_ind;
    int  cxy_nalloc;

    int  *cell;
    int  cell_nalloc;

    int  *a;
    int  *nsubc; /* The number of sub cells for each super cell */
    real *bbcz;  /* Bounding boxes in z for the super cells     */
    real *bb;    /* 3D bounding boxes for the sub cells         */
    int  nc_nalloc;

    int  nsubc_tot;

    gmx_bool print_cycles;
    gmx_nbs_cc_t cc[enbsCCnr];

    gmx_supcell_set_i_x_t *supc_set_i_x;

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
    fprintf(fp,"ns %4d grid %4.1f search %4.1f",
            nbs->cc[enbsCCgrid].count,
            Mcyc_av(&nbs->cc[enbsCCgrid]),
            Mcyc_av(&nbs->cc[enbsCCsearch]));

    if (nbs->nthread_max > 1)
    {
        fprintf(fp," comb %4.1f",
                Mcyc_av(&nbs->cc[enbsCCcombine]));
        fprintf(fp," s. th");
        for(t=0; t<nbs->nthread_max; t++)
        {
            fprintf(fp," %4.1f",
                    Mcyc_av(&nbs->work[t].cc[enbsCCsearch]));
        }
    }
    fprintf(fp,"\n");
}

void gmx_nbsearch_init(gmx_nbsearch_t * nbs_ptr,int natoms_subcell)
{
    gmx_nbsearch_t nbs;
    int t;

    snew(nbs,1);
    *nbs_ptr = nbs;
    
    nbs->naps = natoms_subcell;
    nbs->napc = natoms_subcell*NSUBCELL;

    if (nbs->napc > NBL_NAPC_MAX)
    {
        gmx_fatal(FARGS,
                  "napc (%d) is larger than the maximum allowed value (%d)",
                  nbs->napc,NBL_NAPC_MAX);
    }
    nbs->naps2log = 0;
    while ((1<<nbs->naps2log) < nbs->naps)
    {
        nbs->naps2log++;
    }
    if ((1<<nbs->naps2log) != nbs->naps)
    {
        gmx_fatal(FARGS,"nsbox naps (%d) is not a power of 2",nbs->naps);
    }

    nbs->cxy_na      = NULL;
    nbs->cxy_ind     = NULL;
    nbs->cxy_nalloc  = 0;
    nbs->cell        = NULL;
    nbs->cell_nalloc = 0;
    nbs->a           = NULL;
    nbs->bb          = NULL;
    nbs->nc_nalloc   = 0;

    if (getenv("GMX_NSBOX_BB") != NULL)
    {
        /* Use only bounding box sub cell pair distances,
         * fast, but produces slightly more sub cell pairs.
         */
        nbs->supc_set_i_x = supcell_set_i_x;
        nbs->subc_dc = NULL;
    }
    else
    {
#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
        if (natoms_subcell == 8 && getenv("GMX_NSBOX_NOSSE") == NULL)
        {
            nbs->supc_set_i_x = supcell_set_i_x_sse8;
            nbs->subc_dc = subc_in_range_sse8;
        }
        else
#endif
        {
            nbs->supc_set_i_x = supcell_set_i_x;
            nbs->subc_dc = subc_in_range_x;
        }
    }

    nbs->nthread_max = 1;
#ifdef GMX_OPENMP
    {
        nbs->nthread_max = omp_get_max_threads();
    }
#endif

    /* Initialize the work data structures for each thread */
    snew(nbs->work,nbs->nthread_max);
    for(t=0; t<nbs->nthread_max; t++)
    {
        nbs->work[t].sort_work   = NULL;
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

static int set_grid_size_xy(gmx_nbsearch_t nbs,int n,matrix box)
{
    real adens,tlen,tlen_x,tlen_y,nc_max;
    int  t;

    nbs->atom_density = n/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);

    if (n > nbs->napc)
    {
        /* target cell length */
#if 0
        /* Approximately cubic super cells */
        tlen   = pow(nbs->napc/nbs->atom_density,1.0/3.0);
        tlen_x = tlen;
        tlen_y = tlen;
#else
        /* Approximately cubic sub cells */
        tlen   = pow(nbs->naps/nbs->atom_density,1.0/3.0);
        tlen_x = tlen*NSUBCELL_X;
        tlen_y = tlen*NSUBCELL_Y;
#endif
        /* We round ncx and ncy down, because we get less cell pairs
         * in the nbsist when the fixed cell dimensions (x,y) are
         * larger than the variable one (z) than the other way around.
         */
        nbs->ncx = max(1,(int)(box[XX][XX]/tlen_x));
        nbs->ncy = max(1,(int)(box[YY][YY]/tlen_y));
    }
    else
    {
        nbs->ncx = 1;
        nbs->ncy = 1;
    }

    if (nbs->ncx*nbs->ncy+1 > nbs->cxy_nalloc)
    {
        nbs->cxy_nalloc = over_alloc_large(nbs->ncx*nbs->ncy);
        srenew(nbs->cxy_na,nbs->cxy_nalloc);
        srenew(nbs->cxy_ind,nbs->cxy_nalloc+1);

        for(t=0; t<nbs->nthread_max; t++)
        {
            srenew(nbs->work[t].cxy_na,nbs->cxy_nalloc);
        }
    }

    /* Worst case scenario of 1 atom in each last cell */
    nc_max = n/nbs->napc + nbs->ncx*nbs->ncy;
    if (nc_max > nbs->nc_nalloc)
    {
        int bb_nalloc;

        nbs->nc_nalloc = over_alloc_large(nc_max);
        srenew(nbs->a,nbs->nc_nalloc*nbs->napc);
        srenew(nbs->nsubc,nbs->nc_nalloc);
        srenew(nbs->bbcz,nbs->nc_nalloc*NNBSBB_D);
#ifdef GMX_NBS_BBXXXX
        if (NSUBCELL % SIMD_WIDTH != 0)
        {
            gmx_incons("NSUBCELL is not a multiple of SIMD_WITH");
        }
        bb_nalloc = nbs->nc_nalloc*NSUBCELL/SIMD_WIDTH*NNBSBB_XXXX;
#else  
        bb_nalloc = nbs->nc_nalloc*NSUBCELL*NNBSBB_B;
#endif
        sfree_aligned(nbs->bb);
        /* This snew also zeros the contents, this avoid possible
         * floating exceptions in SSE with the unused bb elements.
         */
        snew_aligned(nbs->bb,bb_nalloc,16);
    }

    nbs->sx = box[XX][XX]/nbs->ncx;
    nbs->sy = box[YY][YY]/nbs->ncy;
    nbs->inv_sx = 1/nbs->sx;
    nbs->inv_sy = 1/nbs->sy;

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

static void print_bbsizes(FILE *fp,const gmx_nbsearch_t nbs)
{
    int  ns,c,s,cs,d;
    dvec ba;

    clear_dvec(ba);
    ns = 0;
    for(c=0; c<nbs->nc; c++)
    {
        for(s=0; s<nbs->nsubc[c]; s++)
        {
            cs = c*NSUBCELL + s;
            for(d=0; d<DIM; d++)
            {
                ba[d] += nbs->bb[cs*NNBSBB_B+NNBSBB_C+d] - nbs->bb[cs*NNBSBB_B+d];
            }
            ns++;
        }
    }
    dsvmul(1.0/ns,ba,ba);

    fprintf(fp,"ns bb: %4.2f %4.2f %4.2f  %4.2f %4.2f %4.2f rel %4.2f %4.2f %4.2f\n",
            nbs->box[XX][XX]/(nbs->ncx*NSUBCELL_X),
            nbs->box[YY][YY]/(nbs->ncy*NSUBCELL_Y),
            nbs->box[ZZ][ZZ]*nbs->ncx*nbs->ncy/(nbs->nc*NSUBCELL_Z),
            ba[XX],ba[YY],ba[ZZ],
            ba[XX]*nbs->ncx*NSUBCELL_X/nbs->box[XX][XX],
            ba[YY]*nbs->ncy*NSUBCELL_Y/nbs->box[YY][YY],
            ba[ZZ]*nbs->nc*NSUBCELL_Z/(nbs->ncx*nbs->ncy*nbs->box[ZZ][ZZ]));
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
    int i,j;

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
        /* Complete the partially filled last cell with copies of the last element.
         * This simplifies the bounding box calculation and avoid
         * numerical issues with atoms that are coincidentally close.
         */
        for(; i<na_round; i++)
        {
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cx);
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cy);
            xnb[j++] = -NBAT_FAR_AWAY*(1 + cz + i);
            j++;
        }
        break;
    default:
        gmx_incons("Unsupported stride");
    }
}
static void sort_columns(gmx_nbsearch_t nbs,
                         int n,rvec *x,
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
    real *xnb;
    real *bb_ptr;

    subdiv_x = nbs->naps;
    subdiv_y = NSUBCELL_X*subdiv_x;
    subdiv_z = NSUBCELL_Y*subdiv_y;

    /* Sort the atoms within each x,y column in 3 dimensions */
    for(cxy=cxy_start; cxy<cxy_end; cxy++)
    {
        cx = cxy/nbs->ncy;
        cy = cxy - cx*nbs->ncy;

        na  = nbs->cxy_na[cxy];
        ncz = nbs->cxy_ind[cxy+1] - nbs->cxy_ind[cxy];
        ash = nbs->cxy_ind[cxy]*nbs->napc;

        /* Sort the atoms within each x,y column on z coordinate */
        sort_atoms(ZZ,FALSE,
                   nbs->a+ash,na,x,
                   0,
                   ncz*nbs->napc*SORT_GRID_OVERSIZE/nbs->box[ZZ][ZZ],
                   ncz*nbs->napc*SGSF,sort_work);

        /* This loop goes over the supercells and subcells along z at once */
        for(sub_z=0; sub_z<ncz*NSUBCELL_Z; sub_z++)
        {
            ash_z = ash + sub_z*subdiv_z;
            na_z  = min(subdiv_z,na-(ash_z-ash));

            /* We have already sorted on z */

            if (sub_z % NSUBCELL_Z == 0)
            {
                cz = sub_z/NSUBCELL_Z;
                c  = nbs->cxy_ind[cxy] + cz ;

                /* The number of atoms in this supercell */
                na_c = min(nbs->napc,na-(ash_z-ash));

                nbs->nsubc[c] = min(NSUBCELL,(na_c+nbs->naps-1)/nbs->naps);

                /* Store the z-boundaries of the super cell */
                nbs->bbcz[c*NNBSBB_D  ] = x[nbs->a[ash_z]][ZZ];
                nbs->bbcz[c*NNBSBB_D+1] = x[nbs->a[ash_z+na_c-1]][ZZ];
            }

#if NSUBCELL_Y > 1
            sort_atoms(YY,(sub_z & 1),
                       nbs->a+ash_z,na_z,x,
                       cy*nbs->sy,nbs->inv_sy,subdiv_y*SGSF,sort_work);
#endif

            for(sub_y=0; sub_y<NSUBCELL_Y; sub_y++)
            {
                ash_y = ash_z + sub_y*subdiv_y;
                na_y  = min(subdiv_y,na-(ash_y-ash));

#if NSUBCELL_X > 1
                sort_atoms(XX,((cz*NSUBCELL_Y + sub_y) & 1),
                           nbs->a+ash_y,na_y,x,
                           cx*nbs->sx,nbs->inv_sx,subdiv_x*SGSF,sort_work);
#endif

                for(sub_x=0; sub_x<NSUBCELL_X; sub_x++)
                {
                    ash_x = ash_y + sub_x*subdiv_x;
                    na_x  = min(subdiv_x,na-(ash_x-ash));

                    xnb   = nbat->x + ash_x*nbat->xstride;

                    if (na_x > 0)
                    {
                        /* Now we have sorted the atoms, set the cell indices */
                        for(a=0; a<na_x; a++)
                        {
                            nbs->cell[nbs->a[ash_x+a]] = ash_x + a;
                        }

                        copy_rvec_to_nbat_real(nbs->a+ash_x,na_x,nbs->naps,x,
                                               nbat->XFormat,xnb,
                                               nbs->naps*(cx*NSUBCELL_X+sub_x),
                                               nbs->naps*(cy*NSUBCELL_Y+sub_y),
                                               nbs->naps*sub_z);

#ifdef GMX_NBS_BBXXXX
                        /* Store the bounding boxes in a format convenient
                         * for SSE calculations: xxxxyyyyzzzz...
                         */
                        bb_ptr =
                            nbs->bb +
                            (ash_x>>(nbs->naps2log+SIMD_WIDTH_2LOG))*NNBSBB_XXXX +
                            ((ash_x>>nbs->naps2log) & (SIMD_WIDTH-1));

                        /* There is something wrong with this
                         * calc_bounding_box_xxxx_sse function or call.
                         * the results are incorrect.
                        if (nbat->XFormat == nbatXYZQ)
                        {
                            calc_bounding_box_xxxx_sse(na_x,xnb,
                                                       nbs->work->bb_tmp,
                                                       bb_ptr);
                        }
                        else
                        */
                        {
                            calc_bounding_box_xxxx(na_x,nbat->xstride,xnb,
                                                   bb_ptr);
                        }
#else
                        /* Store the bounding boxes as xyz.xyz. */
                        bb_ptr = nbs->bb+(ash_x>>nbs->naps2log)*NNBSBB_B;

                        calc_bounding_box(na_x,nbat->xstride,xnb,
                                          bb_ptr);

                        if (gmx_debug_at)
                        {
                            fprintf(debug,"%2d %2d %2d %d %d %d bb %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                                    cx,cy,cz,sub_x,sub_y,sub_z,
                                    (nbs->bb+(ash_x/nbs->naps)*NNBSBB_B)[0],
                                    (nbs->bb+(ash_x/nbs->naps)*NNBSBB_B)[4],
                                    (nbs->bb+(ash_x/nbs->naps)*NNBSBB_B)[1],
                                    (nbs->bb+(ash_x/nbs->naps)*NNBSBB_B)[5],
                                    (nbs->bb+(ash_x/nbs->naps)*NNBSBB_B)[2],
                                    (nbs->bb+(ash_x/nbs->naps)*NNBSBB_B)[6]);
                        }
#endif
                    }
                    else
                    {
                        clear_nbat_real(nbs->naps,nbat->xstride,xnb);
                    }
                }
            }
        }

        /* Set the unused atom indices to -1 */
        for(ind=na; ind<ncz*nbs->napc; ind++)
        {
            nbs->a[ash+ind] = -1;
        }

        /*
        copy_rvec_to_nbat_real(axy,na,ncz*nbs->napc,x,nbat->xstride,xnb);

        calc_bounding_box(ncz,nbs->napc,nbat->xstride,xnb,
                          nbs->bb+nbs->cxy_ind[i]*NNBSBB);
        */
    }
}

static void calc_cell_indices(gmx_nbsearch_t nbs,
                              int n,rvec *x,
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

        for(i=0; i<nbs->ncx*nbs->ncy; i++)
        {
            cxy_na[i] = 0;
        }

        n0 = (int)((t+0)*n)/nthread;
        n1 = (int)((t+1)*n)/nthread;
        for(i=n0; i<n1; i++)
        {
            /* We need to be careful with rounding,
             * particles might be a few bits outside the local box.
             * The int cast takes care of the lower bound,
             * we need to explicitly take care of the upper bound.
             */
            cx = (int)(x[i][XX]*nbs->inv_sx);
            if (cx == nbs->ncx)
            {
                cx = nbs->ncx - 1;
            }
            cy = (int)(x[i][YY]*nbs->inv_sy);
            if (cy == nbs->ncy)
            {
                cy = nbs->ncy - 1;
            }
            /* For the moment cell contains only the x and y indices, not z */
            nbs->cell[i] = cx*nbs->ncy + cy;
            cxy_na[nbs->cell[i]]++;
        }
    }

    /* Make the cell index as a function of x and y */
    ncz_max = 0;
    nbs->cxy_ind[0] = 0;
    for(i=0; i<nbs->ncx*nbs->ncy; i++)
    {
        cxy_na_i = nbs->work[0].cxy_na[i];
        for(t=1; t<nthread; t++)
        {
            cxy_na_i += nbs->work[t].cxy_na[i];
        }
        ncz = (cxy_na_i + nbs->napc - 1)/nbs->napc;
        nbs->cxy_ind[i+1] = nbs->cxy_ind[i] + ncz;
        if (ncz > ncz_max)
        {
            ncz_max = ncz;
        }
        /* Clear cxy_na, so we can reuse the array below */
        nbs->cxy_na[i] = 0;
    }
    nbs->nc = nbs->cxy_ind[nbs->ncx*nbs->ncy];

    nbat->natoms = nbs->nc*nbs->napc;

    if (debug)
    {
        fprintf(debug,"ns napc %d naps %d super-cells: %d x %d y %d z %.1f maxz %d\n",
                nbs->napc,nbs->naps,nbs->nc,
                nbs->ncx,nbs->ncy,nbs->nc/((double)(nbs->ncx*nbs->ncy)),
                ncz_max);
        if (gmx_debug_at)
        {
            i = 0;
            for(cy=0; cy<nbs->ncy; cy++)
            {
                for(cx=0; cx<nbs->ncx; cx++)
                {
                    fprintf(debug," %2d",nbs->cxy_ind[i+1]-nbs->cxy_ind[i]);
                    i++;
                }
                fprintf(debug,"\n");
            }
        }
    }

    /* Make sure the work array for sorting is large enough */
    if (ncz_max*nbs->napc*SGSF > nbs->work[0].sort_work_nalloc)
    {
        for(t=0; t<nbs->nthread_max; t++)
        {
            nbs->work[t].sort_work_nalloc =
                over_alloc_large(ncz_max*nbs->napc*SGSF);
            srenew(nbs->work[t].sort_work,nbs->work[t].sort_work_nalloc);
        }
    }

    /* Now we know the dimensions we can fill the grid.
     * This is the first, unsorted fill. We sort the columns after this.
     */
    for(i=0; i<n; i++)
    {
        cxy = nbs->cell[i];
        nbs->a[nbs->cxy_ind[cxy]*nbs->napc + nbs->cxy_na[cxy]++] = i;
    }

    /* Sort the super-cell columns along z into the sub-cells. */
    nthread = omp_get_max_threads();
#pragma omp parallel for schedule(static)
    for(t=0; t<nthread; t++)
    {
        sort_columns(nbs,n,x,nbat,
                     ((t+0)*nbs->ncx*nbs->ncy)/nthread,
                     ((t+1)*nbs->ncx*nbs->ncy)/nthread,
                     nbs->work[t].sort_work);
    }

    if (debug)
    {
        int c;

        nbs->nsubc_tot = 0;
        for(c=0; c<nbs->nc; c++)
        {
            nbs->nsubc_tot += nbs->nsubc[c];
        }
        fprintf(debug,"ns non-zero sub-cells: %d average atoms %.2f\n",
                nbs->nsubc_tot,n/(double)nbs->nsubc_tot);

        print_bbsizes(debug,nbs);
    }
}

static void nb_realloc_void(void **ptr,
                            int nbytes_copy,int nbytes_new,
                            gmx_nbat_alloc_t *ma,
                            gmx_nbat_free_t  *mf)
{
    void *ptr_new;

    ma(&ptr_new,nbytes_new);

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
    nb_realloc_int(&nbat->type,n,nbat->alloc,nbat->free);
    if (nbat->XFormat != nbatXYZQ)
    {
        nb_realloc_real(&nbat->q,n,nbat->alloc,nbat->free);
    }
    nb_realloc_real(&nbat->x,n*nbat->xstride,nbat->alloc,nbat->free);
    nb_realloc_real(&nbat->f,n*nbat->xstride,nbat->alloc,nbat->free);
    nbat->nalloc = n;
}

void gmx_nbsearch_put_on_grid(gmx_nbsearch_t nbs,
                              int ePBC,matrix box,int n,rvec *x,
                              gmx_nb_atomdata_t *nbat)
{
    int nc_max;

    nbs_cycle_start(&nbs->cc[enbsCCgrid]);

    nbs->ePBC = ePBC;
    copy_mat(box,nbs->box);

    nc_max = set_grid_size_xy(nbs,n,nbs->box);

    if (n > nbs->cell_nalloc)
    {
        nbs->cell_nalloc = over_alloc_large(n);
        srenew(nbs->cell,nbs->cell_nalloc);
    }

    if (nc_max*nbs->napc > nbat->nalloc)
    {
        gmx_nb_atomdata_realloc(nbat,nc_max*nbs->napc);
    }

    calc_cell_indices(nbs,n,x,nbat);

    nbs_cycle_stop(&nbs->cc[enbsCCgrid]);
}

static void get_cell_range(real b0,real b1,int nc,real s,real invs,
                           real d2,real r2,int *cf,int *cl)
{
    *cf = max((int)(b0*invs),0);
    
    while (*cf > 0 && d2 + sqr(b0 - (*cf-1+1)*s) < r2)
    {
        (*cf)--;
    }

    *cl = min((int)(b1*invs),nc-1);
    while (*cl < nc-1 && d2 + sqr((*cl+1)*s - b1) < r2)
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

static void check_subcell_list_space(gmx_nblist_t *nbl,int nsupercell)
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
    *ptr = save_calloc_aligned("ptr",__FILE__,__LINE__,nbytes,1,16);
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

static void print_nblist_statistics(FILE *fp,const gmx_nblist_t *nbl,
                                    const gmx_nbsearch_t nbs,real rl)
{
    int i,j4,j,si,b;
    int c[NSUBCELL+1];

    fprintf(fp,"nbl nci %d nsj4 %d nsi %d excl4 %d\n",
            nbl->nci,nbl->nsj4,nbl->nsi,nbl->nexcl);
    fprintf(fp,"nbl naps %d rl %g ncp %d per cell %.1f atoms %.1f ratio %.2f\n",
            nbl->naps,rl,nbl->nsi,nbl->nsi/(double)nbs->nsubc_tot,
            nbl->nsi/(double)nbs->nsubc_tot*nbs->naps,
            nbl->nsi/(double)nbs->nsubc_tot*nbs->naps/(0.5*4.0/3.0*M_PI*rl*rl*rl*nbs->nsubc_tot*nbs->naps/det(nbs->box)));

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

static void make_subcell_list(const gmx_nbsearch_t nbs,
                              gmx_nblist_t *nbl,
                              int ci,int cj,
                              gmx_bool ci_equals_cj,
                              int stride,const real *x,
                              real rl2,real rbb2)
{
    int  naps;
    int  npair;
    int  sj,si1,si,csj,csi;
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

    naps = nbs->naps;

    for(sj=0; sj<nbs->nsubc[cj]; sj++)
    {
        sj4_ind   = (nbl->work->sj_ind >> 2);
        sj_offset = nbl->work->sj_ind - sj4_ind*4;
        sj4       = &nbl->sj4[sj4_ind];
        
        csj = cj*NSUBCELL + sj;

        /* Initialize this j-subcell i-subcell list */
        sj4->sj[sj_offset] = csj;
        imask              = 0;

        if (!nbl->TwoWay && ci_equals_cj)
        {
            si1 = sj + 1;
        }
        else
        {
            si1 = nbs->nsubc[ci];
        }

#ifdef GMX_NBS_BBXXXX
        /* Determine all si1 bb distances in one call with SSE */
        subc_bb_dist2_sse_xxxx(nbs->bb+(csj>>SIMD_WIDTH_2LOG)*NNBSBB_XXXX+(csj & (SIMD_WIDTH-1)),
                               si1,bb_ci,d2l);
#endif

        npair = 0;
        for(si=0; si<si1; si++)
        {
#ifndef GMX_NBS_BBXXXX
            /* Determine the bb distance between csi and csj */
            d2l[si] = subc_bb_dist2(naps,si,bb_ci,csj,nbs->bb);
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
                (d2 < rl2 && nbs->subc_dc(naps,si,x_ci,csj,stride,x,rl2)))
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
            if (!nbs->subc_dc(naps,si_last,x_ci,csj,stride,x,rl2))
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
                         const gmx_nbl_ci_t *nbl_ci,
                         const t_blocka *excl)
{
    const int *cell;
    int naps;
    int ci;
    int sj_ind_last;
    int sj_first,sj_last;
    int ndirect;
    int i,ai,si,eind,ge,se;
    int found,sj_ind_0,sj_ind_1,sj_ind_m;
    int sj_m;
    gmx_bool Found_si;
    int si_ind;
    gmx_nbl_excl_t *nbl_excl;
    int inner_i,inner_e,w;

    cell = nbs->cell;

    naps = nbs->naps;

    if (nbl_ci->sj4_ind_end == nbl_ci->sj4_ind_start)
    {
        /* Empty list */
        return;
    }

    ci = nbl_ci->ci;

    sj_ind_last = nbl->work->sj_ind - 1;

    sj_first = nbl->sj4[nbl_ci->sj4_ind_start].sj[0];
    sj_last  = nbl_sj(nbl,sj_ind_last);

    /* Determine how many contiguous j-sub-cells we have starting
     * from the first i-sub-cell. This number can be used to directly
     * calculate j-sub-cell indices for excluded atoms.
     */
    ndirect = 0;
    while (nbl_ci->sj4_ind_start*4 + ndirect <= sj_ind_last &&
           nbl_sj(nbl,nbl_ci->sj4_ind_start*4+ndirect) == ci*NSUBCELL + ndirect)
    {
        ndirect++;
    }

    /* Loop over the atoms in the i super-cell */
    for(i=0; i<nbs->napc; i++)
    {
        ai = nbs->a[ci*nbs->napc+i];
        if (ai >= 0)
        {
            si  = (i>>nbs->naps2log);

#ifdef DEBUG_NSBOX_EXCLS
            if (nbs->cell[ai] != ci*nbs->napc + i)
            {
                gmx_incons("Index mismatch");
            }
#endif

            /* Loop over the exclusions for this i-atom */
            for(eind=excl->index[ai]; eind<excl->index[ai+1]; eind++)
            {
                ge = cell[excl->a[eind]];

                /* Without shifts we only calculate interactions j>i
                 * for one-way neighbor lists.
                 */
                if (!nbl->TwoWay && nbl_ci->shift == CENTRAL &&
                    ge <= ci*nbs->napc + i)
                {
                    continue;
                }

                se = ge>>nbs->naps2log;
                /* Could the sub-cell se be in our list? */
                if (se >= sj_first && se <= sj_last)
                {
                    if (se < sj_first + ndirect)
                    {
                        /* We can calculate sj_ind directly from se */
                        found = nbl_ci->sj4_ind_start*4 + se - sj_first;
                    }
                    else
                    {
                        /* Search for se using bisection */
                        found = -1;
                        sj_ind_0 = nbl_ci->sj4_ind_start*4 + ndirect;
                        sj_ind_1 = sj_ind_last + 1;
                        while (found == -1 && sj_ind_0 < sj_ind_1)
                        {
                            sj_ind_m = (sj_ind_0 + sj_ind_1)>>1;
                            sj_m = nbl_sj(nbl,sj_ind_m);
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

                    if (found >= 0 &&
                        (nbl_imask0(nbl,found) & (1U << ((found & 3)*NSUBCELL + si))))
                    {
                        inner_i = i  - si*naps;
                        inner_e = ge - se*naps;
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

static void nb_realloc_ci(gmx_nblist_t *nbl,int n)
{
    nbl->ci_nalloc = over_alloc_small(n);
    nb_realloc_void((void **)&nbl->ci,
                    nbl->nci*sizeof(*nbl->ci),
                    nbl->ci_nalloc*sizeof(*nbl->ci),
                    nbl->alloc,nbl->free);
}

static void new_ci_entry(gmx_nblist_t *nbl,int ci,int shift,
                         gmx_nbl_work_t *work)
{
    if (nbl->nci + 1 > nbl->ci_nalloc)
    {
        nb_realloc_ci(nbl,nbl->nci+1);
    }
    nbl->ci[nbl->nci].ci            = ci;
    nbl->ci[nbl->nci].shift         = shift;
    nbl->ci[nbl->nci].sj4_ind_start = nbl->nsj4;
    nbl->ci[nbl->nci].sj4_ind_end   = nbl->nsj4;
}

static void close_ci_entry(gmx_nblist_t *nbl,int max_j4list_av,int nc_bal)
{
    int j4len,tlen;
    int nb,b;
    int max_j4list;

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

        if (max_j4list_av > 0)
        {
            /* The first ci blocks should be larger, to avoid overhead.
             * The last ci blocks should be smaller, to improve load balancing.
             */
            max_j4list = max(1,
                             max_j4list_av*nc_bal*3/(2*(nbl->nci - 1 + nc_bal)));

            if (j4len > max_j4list)
            {
                /* Split ci in the minimum number of blocks <=jlen */
                nb = (j4len + max_j4list - 1)/max_j4list;
                /* Make blocks similar sized, last one smallest */
                tlen = (j4len + nb - 1)/nb;
                
                if (nbl->nci + nb - 1 > nbl->ci_nalloc)
                {
                    nb_realloc_ci(nbl,nbl->nci+nb-1);
                }
                
                /* Set the end of the last block to the current end */
                nbl->ci[nbl->nci+nb-2].sj4_ind_end =
                    nbl->ci[nbl->nci-1].sj4_ind_end;

                for(b=1; b<nb; b++)
                {
                    nbl->ci[nbl->nci-1].sj4_ind_end =
                        nbl->ci[nbl->nci-1].sj4_ind_start + tlen;
                    nbl->ci[nbl->nci].sj4_ind_start =
                        nbl->ci[nbl->nci-1].sj4_ind_end;
                    nbl->ci[nbl->nci].ci            = nbl->ci[nbl->nci-1].ci;
                    nbl->ci[nbl->nci].shift         = nbl->ci[nbl->nci-1].shift;
                    nbl->nci++;
                }
            }
        }
    }
}

static void clear_nblist(gmx_nblist_t *nbl)
{
    nbl->nci          = 0;
    nbl->nsj4         = 0;
    nbl->work->sj_ind = nbl->nsj4*4;
    nbl->nsi          = 0;
    nbl->nexcl        = 1;

    nbl->work->sj4_init = 0;
}

static void set_supcell_i_bb(const real *bb,int ci,
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

static void supcell_set_i_x(const gmx_nbsearch_t nbs,int ci,
                            real shx,real shy,real shz,
                            int stride,const real *x,
                            gmx_nbl_work_t *work)
{
    int  ia,i;
    real *x_ci;

    x_ci = work->x_ci;

    ia = ci*NSUBCELL*nbs->naps;
    for(i=0; i<nbs->napc; i++)
    {
        x_ci[i*3 + 0] = x[(ia+i)*stride + 0] + shx;
        x_ci[i*3 + 1] = x[(ia+i)*stride + 1] + shy;
        x_ci[i*3 + 2] = x[(ia+i)*stride + 2] + shz;
    }
}

static void supcell_set_i_x_sse8(const gmx_nbsearch_t nbs,int ci,
                                 real shx,real shy,real shz,
                                 int stride,const real *x,
                                 gmx_nbl_work_t *work)
{
    int  si,io,ia,i,j;
    real *x_ci;

    x_ci = work->x_ci;

    for(si=0; si<NSUBCELL; si++)
    {
        for(i=0; i<nbs->naps; i+=4)
        {
            io = si*nbs->naps + i;
            ia = ci*NSUBCELL*nbs->naps + io;
            for(j=0; j<4; j++)
            {
                x_ci[io*3 + j + 0] = x[(ia+j)*stride+0] + shx;
                x_ci[io*3 + j + 4] = x[(ia+j)*stride+1] + shy;
                x_ci[io*3 + j + 8] = x[(ia+j)*stride+2] + shz;
            }
        }
    }
}

static int get_max_j4list(const gmx_nbsearch_t nbs,
                          gmx_nblist_t *nbl,
                          int min_ci_balanced)
{
    real xy_diag,r_eff_sup;
    int  nj4_est,nparts;
    int  max_j4list;

    /* The average diagonal of a super cell */
    xy_diag = sqrt(sqr(nbs->box[XX][XX]/nbs->ncx) +
                   sqr(nbs->box[YY][YY]/nbs->ncy) +
                   sqr(nbs->box[ZZ][ZZ]*nbs->ncx*nbs->ncy/nbs->nc));

    /* The formulas below are a heuristic estimate of the average nj per ci*/
    r_eff_sup = nbl->rlist + 0.4*xy_diag;
    
    nj4_est = (int)(0.5*4.0/3.0*M_PI*pow(r_eff_sup,3)*nbs->atom_density/(4*nbs->naps) + 0.5);

    if (min_ci_balanced <= 0 || nbs->nc >= min_ci_balanced)
    {
        /* We don't need to worry */
        max_j4list = -1;
    }
    else
    {
        /* Estimate the number of parts we need to cut each full list
         * for one i super cell into.
         */
        nparts = (min_ci_balanced + nbs->nc - 1)/nbs->nc;
        /* Thus the (average) maximum j-list size should be as follows */
        max_j4list = max(1,(nj4_est + nparts - 1)/nparts);
    }

    if (debug)
    {
        fprintf(debug,"nbl nj4 estimate %d, max_j4list %d\n",
                nj4_est,max_j4list);
    }

    return max_j4list;
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
    int n,i,j4,j0,j1;
    int sj4_offset,si_offset,excl_offset;
    const gmx_nblist_t *nbli;
    int nth,th;

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

    nth = omp_get_max_threads();

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
        /* The ci list copy is probably not work parallelizing */
        for(i=0; i<nbli->nci; i++)
        {
            nblc->ci[nblc->nci]                = nbli->ci[i];
            nblc->ci[nblc->nci].sj4_ind_start += sj4_offset;
            nblc->ci[nblc->nci].sj4_ind_end   += sj4_offset;
            nblc->nci++;
        }
#pragma omp parallel private(th,j0,j1,j4)
        {
            th = omp_get_thread_num();
            j0 = (nbli->nsj4*(th+0))/nth;
            j1 = (nbli->nsj4*(th+1))/nth;
            for(j4=j0; j4<j1; j4++)
            {
                nblc->sj4[nblc->nsj4+j4] = nbli->sj4[j4];
                nblc->sj4[nblc->nsj4+j4].imei[0].excl_ind += excl_offset;
                nblc->sj4[nblc->nsj4+j4].imei[1].excl_ind += excl_offset;
            }

            j0 = (nbli->nexcl*(th+0))/nth;
            j1 = (nbli->nexcl*(th+1))/nth;
            for(j4=j0; j4<j1; j4++)
            {
                nblc->excl[nblc->nexcl+j4] = nbli->excl[j4];
            }
            //memcpy(nbli->excl,nblc->excl+nblc->nexcl+j0,(j1-j0)*sizeof(*nblc->excl));
        }

        nblc->nsj4  += nbli->nsj4;
        nblc->nsi   += nbli->nsi;
        nblc->nexcl += nbli->nexcl;
    }
}

static void gmx_nbsearch_make_nblist_part(const gmx_nbsearch_t nbs,
                                          gmx_nbs_work_t *work,
                                          const gmx_nb_atomdata_t *nbat,
                                          const t_blocka *excl,
                                          real rcut,real rlist,
                                          int min_ci_balanced,
                                          int ci_start,int ci_end,
                                          gmx_nblist_t *nbl)
{
    gmx_bool bDomDec;
    int  max_j4list;
    matrix box;
    real rl2,rbb2;
    int  d,ci,ci_xy,ci_x,ci_y,cj;
    ivec shp;
    int  tx,ty,tz;
    int  shift;
    gmx_bool bMakeList;
    real shx,shy,shz;
    real *bbcz,*bb;
    real bx0,bx1,by0,by1,bz0,bz1;
    real bz1_frac;
    real d2z,d2zx,d2zxy,d2xy;
    int  cxf,cxl,cyf,cyf_x,cyl;
    int  cx,cy;
    int  c0,c1,cs,cf,cl;

    nbs_cycle_start(&work->cc[enbsCCsearch]);

    bDomDec = FALSE;

    nbl->napc = nbs->napc;
    nbl->naps = nbs->naps;

    /* Currently this code only makes two-way lists */
    nbl->TwoWay = FALSE;

    nbl->rcut   = rcut;
    nbl->rlist  = rlist;

    max_j4list = get_max_j4list(nbs,nbl,min_ci_balanced);

    clear_nblist(nbl);

    copy_mat(nbs->box,box);

    rl2 = nbl->rlist*nbl->rlist;

    if (nbs->subc_dc == NULL)
    {
        rbb2 = rl2;
    }
    else
    {
        /* If the distance between two sub-cell bounding boxes is less
         * than the nblist cut-off minus half of the average x/y diagonal
         * spacing of the sub-cells, do not check the distance between
         * all particle pairs in the sub-cell, since this pairs is very
         * likely to have atom pairs within the cut-off.
         */
        rbb2 = sqr(max(0,
                       nbl->rlist -
                       0.5*sqrt(sqr(box[XX][XX]/(nbs->ncx*NSUBCELL_X)) +
                                sqr(box[YY][YY]/(nbs->ncy*NSUBCELL_Y)))));
    }
    if (debug)
    {
        fprintf(debug,"nbl bounding box only distance %f\n",sqrt(rbb2));
    }

    /* Set the shift range */
    for(d=0; d<DIM; d++)
    {
        /* We need to add these domain shift limits for DD
        sh0[d] = -1;
        sh1[d] = 1;
        */
        /* Check if we need periodicity shifts.
         * Without PBC or with domain decomposition we don't need them.
         */
        /*
        if (d >= ePBC2npbcdim(fr->ePBC) || (bDomDec && dd->nc[d] > 1))
        */
        if (d >= ePBC2npbcdim(nbs->ePBC))
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

    bbcz = nbs->bbcz;
    bb   = nbs->bb;

    ci_xy = 0;
    for(ci=ci_start; ci<ci_end; ci++)
    {
        while (ci >= nbs->cxy_ind[ci_xy+1])
        {
            ci_xy++;
        }
        ci_x = ci_xy/nbs->ncy;
        ci_y = ci_xy - ci_x*nbs->ncy;

        /* Loop over shift vectors in three dimensions */
        for (tz=-shp[ZZ]; tz<=shp[ZZ]; tz++)
        {
            shz = tz*box[ZZ][ZZ];

            bz0 = bbcz[ci*NNBSBB_D  ] + shz;
            bz1 = bbcz[ci*NNBSBB_D+1] + shz;

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

            if (d2z >= rl2)
            {
                continue;
            }

            bz1_frac =
                bz1/((real)(nbs->cxy_ind[ci_xy+1] - nbs->cxy_ind[ci_xy]));
            if (bz1_frac < 0)
            {
                bz1_frac = 0;
            }
            /* The check with bz1_frac close to or larger than 1 comes later */

            for (ty=-shp[YY]; ty<=shp[YY]; ty++)
            {
                shy = ty*box[YY][YY] + tz*box[ZZ][YY];
            
                by0 = (ci_y  )*nbs->sy + shy;
                by1 = (ci_y+1)*nbs->sy + shy;

                get_cell_range(by0,by1,nbs->ncy,nbs->sy,nbs->inv_sy,d2z,rl2,
                               &cyf,&cyl);

                for (tx=-shp[XX]; tx<=shp[XX]; tx++)
                {
                    shift = XYZ2IS(tx,ty,tz);

#ifdef NSBOX_SHIFT_BACKWARD
                    if (shift > CENTRAL)
                    {
                        continue;
                    }
#endif

                    shx = tx*box[XX][XX] + ty*box[YY][XX] + tz*box[ZZ][XX];

                    bx0 = (ci_x  )*nbs->sx + shx;
                    bx1 = (ci_x+1)*nbs->sx + shx;

                    get_cell_range(bx0,bx1,nbs->ncx,nbs->sx,nbs->inv_sx,d2z,rl2,
                                   &cxf,&cxl); 

                    new_ci_entry(nbl,ci,shift,nbl->work);

#ifndef NSBOX_SHIFT_BACKWARD
                    if (!nbl->TwoWay && cxf < ci_x)
#else
                    if (!nbl->TwoWay && shift == CENTRAL && cxf < ci_x)
#endif
                    {
                        /* Leave the pairs with i > j.
                         * x is the major index, so skip half of it.
                         */
                        cxf = ci_x;
                    }

                    set_supcell_i_bb(nbs->bb,ci,shx,shy,shz,nbl->work->bb_ci);

                    nbs->supc_set_i_x(nbs,ci,shx,shy,shz,
                                      nbat->xstride,nbat->x,nbl->work);

                    for(cx=cxf; cx<=cxl; cx++)
                    {
                        d2zx = d2z;
                        if (cx*nbs->sx > bx1)
                        {
                            d2zx += sqr(cx*nbs->sx - bx1);
                        }
                        else if ((cx+1)*nbs->sx < bx0)
                        {
                            d2zx += sqr((cx+1)*nbs->sx - bx0);
                        }

#ifndef NSBOX_SHIFT_BACKWARD
                        if (!nbl->TwoWay && cx == 0 && cyf < ci_y)
#else
                        if (!nbl->TwoWay &&
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
                            c0 = nbs->cxy_ind[cx*nbs->ncy+cy];
                            c1 = nbs->cxy_ind[cx*nbs->ncy+cy+1];
#ifdef NSBOX_SHIFT_BACKWARD
                            if (!nbl->TwoWay && shift == CENTRAL && c0 < ci)
                            {
                                c0 = ci;
                            }
#endif

                            d2zxy = d2zx;
                            if (cy*nbs->sy > by1)
                            {
                                d2zxy += sqr(cy*nbs->sy - by1);
                            }
                            else if ((cy+1)*nbs->sy < by0)
                            {
                                d2zxy += sqr((cy+1)*nbs->sy - by0);
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
                                      (bbcz[cf*NNBSBB_D+1] >= bz0 ||
                                       d2xy + sqr(bbcz[cf*NNBSBB_D+1] - bz0) < rl2))
                                {
                                    cf--;
                                }

                                /* Find the highest cell that can possibly
                                 * be within range.
                                 */
                                cl = cs;
                                while(cl < c1-1 &&
                                      (bbcz[cl*NNBSBB_D] <= bz1 ||
                                       d2xy + sqr(bbcz[cl*NNBSBB_D] - bz1) < rl2))
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

                                if (cf <= cl)
                                {
                                    check_subcell_list_space(nbl,cl-cf+1);

                                    for(cj=cf; cj<=cl; cj++)
                                    {
                                        make_subcell_list(nbs,nbl,ci,cj,
                                                          (shift == CENTRAL && ci == cj),
                                                          nbat->xstride,nbat->x,
                                                          rl2,rbb2);
                                    }
                                }
                            }
                        }  
                    }

                    /* Set the exclusions for this ci list */
                    set_ci_excls(nbs,nbl,&(nbl->ci[nbl->nci]),excl);

                    /* Close this ci list */
                    close_ci_entry(nbl,max_j4list,min_ci_balanced);
                }
            }
        }
    }

    nbs_cycle_stop(&work->cc[enbsCCsearch]);
    
    if (debug)
    {
        print_nblist_statistics(debug,nbl,nbs,rlist);

        if (gmx_debug_at)
        {
            print_nblist_ci_sj(debug,nbl);
        }
    }
}

void gmx_nbsearch_make_nblist(const gmx_nbsearch_t nbs,
                              const gmx_nb_atomdata_t *nbat,
                              const t_blocka *excl,
                              real rcut,real rlist,
                              int min_ci_balanced,
                              int nnbl,gmx_nblist_t **nbl,
                              gmx_bool CombineNBLists)
{
    int nth,th;

    if (debug)
    {
        fprintf(debug,"ns making %d nblists\n",nnbl);
    }
    nbs_cycle_start(&nbs->cc[enbsCCsearch]);

#pragma omp parallel for schedule(static)
    for(th=0; th<nnbl; th++)
    {
        /* Divide the i super cell equally over the nblists */
        gmx_nbsearch_make_nblist_part(nbs,&nbs->work[th],nbat,excl,
                                      rcut,rlist,min_ci_balanced,
                                      ((th+0)*nbs->nc)/nnbl,
                                      ((th+1)*nbs->nc)/nnbl,
                                      nbl[th]);
    }
    nbs_cycle_stop(&nbs->cc[enbsCCsearch]);

    if (CombineNBLists && nnbl > 1)
    {
        nbs_cycle_start(&nbs->cc[enbsCCcombine]);

        combine_nblists(nnbl-1,nbl+1,nbl[0]);
 
        nbs_cycle_stop(&nbs->cc[enbsCCcombine]);
    }

    if (nbs->print_cycles &&
        nbs->cc[enbsCCgrid].count > 0 &&
        nbs->cc[enbsCCgrid].count % 10 == 0)
    {
        nbs_cycle_print(stderr,nbs);
    }

    if (debug)
    {
        print_nblist_statistics(debug,nbl[0],nbs,rlist);
    }
    if (gmx_debug_at)
    {
        print_nblist_ci_sj(debug,nbl[0]);
    }
}

void gmx_nb_atomdata_init(gmx_nb_atomdata_t *nbat,int nbatXFormat,
                          int ntype,const real *nbfp,
                          gmx_nbat_alloc_t *alloc,
                          gmx_nbat_free_t  *free)
{
    int i,j;

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
    for(i=0; i<nbat->ntype; i++)
    {
        for(j=0; j<nbat->ntype; j++)
        {
            if (i < ntype && j < ntype)
            {
                nbat->nbfp[(i*nbat->ntype+j)*2  ] = nbfp[(i*ntype+j)*2  ];
                nbat->nbfp[(i*nbat->ntype+j)*2+1] = nbfp[(i*ntype+j)*2+1];
            }
            else
            {
                /* Add zero parameters for the additional dummy atom type */
                nbat->nbfp[(i*nbat->ntype+j)*2  ] = 0;
                nbat->nbfp[(i*nbat->ntype+j)*2+1] = 0;
            }
        }
    }

    nbat->natoms  = 0;
    nbat->type    = NULL;
    nbat->XFormat = nbatXFormat;
    nbat->q       = NULL;
    nbat->alloc((void **)&nbat->shift_vec,SHIFTS*sizeof(*nbat->shift_vec));
    nbat->xstride = (nbatXFormat == nbatXYZQ ? 4 : 3);
    nbat->x       = NULL;
    nbat->nalloc  = 0;
}

void gmx_nb_atomdata_set_atomtypes(gmx_nb_atomdata_t *nbat,
                                   const gmx_nbsearch_t nbs,
                                   const int *type)
{
    int i,ncz,ash;

    /* Loop over all columns and copy and fill */
    for(i=0; i<nbs->ncx*nbs->ncy; i++)
    {
        ncz = nbs->cxy_ind[i+1] - nbs->cxy_ind[i];
        ash = nbs->cxy_ind[i]*nbs->napc;

        copy_int_to_nbat_int(nbs->a+ash,nbs->cxy_na[i],ncz*nbs->napc,
                             type,nbat->ntype-1,nbat->type+ash);
    }
}

void gmx_nb_atomdata_set_charges(gmx_nb_atomdata_t *nbat,
                                 const gmx_nbsearch_t nbs,
                                 const real *charge)
{
    int  cxy,ncz,ash,na,na_round,i,j;
    real *q;

    /* Loop over all columns and copy and fill */
    for(cxy=0; cxy<nbs->ncx*nbs->ncy; cxy++)
    {
        ash = nbs->cxy_ind[cxy]*nbs->napc;
        na  = nbs->cxy_na[cxy];
        na_round = (nbs->cxy_ind[cxy+1] - nbs->cxy_ind[cxy])*nbs->napc;

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
                                      rvec *x,
                                      gmx_nb_atomdata_t *nbat)
{
    int cxy,na,ash;

#pragma omp parallel for schedule(static) private(na,ash)
    for(cxy=0; cxy<nbs->ncx*nbs->ncy; cxy++)
    {
        na  = nbs->cxy_na[cxy];
        ash = nbs->cxy_ind[cxy]*nbs->napc;

        /* We fill only the real particle locations.
         * We assume the filling entries at the end have been
         * properly set before during ns.
         */
        copy_rvec_to_nbat_real(nbs->a+ash,na,na,x,
                               nbat->XFormat,nbat->x+ash*nbat->xstride,
                               0,0,0);
    }
}

static void gmx_nb_atomdata_add_nbat_f_to_f_part(const gmx_nbsearch_t nbs,
                                                 const gmx_nb_atomdata_t *nbat,
                                                 int a0,int a1,
                                                 rvec *f)
{
    int  a,i;
    const int  *cell;
    const real *fnb;

    cell = nbs->cell;

    fnb = nbat->f;

    /* Loop over all columns and copy and fill */
    for(a=a0; a<a1; a++)
    {
        i = cell[a]*nbat->xstride;

        f[a][XX] += fnb[i];
        f[a][YY] += fnb[i+1];
        f[a][ZZ] += fnb[i+2];
    }
}

void gmx_nb_atomdata_add_nbat_f_to_f(const gmx_nbsearch_t nbs,
                                     const gmx_nb_atomdata_t *nbat,
                                     int natoms,rvec *f)
{
    int nth,th;

    nth = omp_get_max_threads();
#pragma omp parallel for schedule(static)
    for(th=0; th<nth; th++)
    {
        gmx_nb_atomdata_add_nbat_f_to_f_part(nbs,nbat,
                                             ((th+0)*natoms)/nth,
                                             ((th+1)*natoms)/nth,
                                             f);
    }
}
