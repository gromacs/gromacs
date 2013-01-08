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

#ifndef _nbnxn_internal_h
#define _nbnxn_internal_h

#include "typedefs.h"
#include "domdec.h"
#include "gmx_cyclecounter.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef GMX_X86_SSE2
/* Use 4-way SIMD for, always, single precision bounding box calculations */
#define NBNXN_SEARCH_BB_SSE
#endif


#ifdef GMX_NBNXN_SIMD
/* Memory alignment in bytes as required by SIMD aligned loads/stores */
#define NBNXN_MEM_ALIGN  (GMX_NBNXN_SIMD_BITWIDTH/8)
#else
/* No alignment required, but set it so we can call the same routines */
#define NBNXN_MEM_ALIGN  32
#endif


/* A pair-search grid struct for one domain decomposition zone */
typedef struct {
    rvec c0;             /* The lower corner of the (local) grid        */
    rvec c1;             /* The upper corner of the (local) grid        */
    real atom_density;   /* The atom number density for the local grid  */

    gmx_bool bSimple;    /* Is this grid simple or super/sub            */
    int  na_c;           /* Number of atoms per cluster                 */
    int  na_cj;          /* Number of atoms for list j-clusters         */
    int  na_sc;          /* Number of atoms per super-cluster           */
    int  na_c_2log;      /* 2log of na_c                                */

    int  ncx;            /* Number of (super-)cells along x             */
    int  ncy;            /* Number of (super-)cells along y             */
    int  nc;             /* Total number of (super-)cells               */

    real sx;             /* x-size of a (super-)cell                    */
    real sy;             /* y-size of a (super-)cell                    */
    real inv_sx;         /* 1/sx                                        */
    real inv_sy;         /* 1/sy                                        */

    int  cell0;          /* Index in nbs->cell corresponding to cell 0  */

    int  *cxy_na;        /* The number of atoms for each column in x,y  */
    int  *cxy_ind;       /* Grid (super)cell index, offset from cell0   */
    int  cxy_nalloc;     /* Allocation size for cxy_na and cxy_ind      */

    int   *nsubc;        /* The number of sub cells for each super cell */
    float *bbcz;         /* Bounding boxes in z for the super cells     */
    float *bb;           /* 3D bounding boxes for the sub cells         */
    float *bbj;          /* 3D j-b.boxes for SSE-double or AVX-single   */
    int   *flags;        /* Flag for the super cells                    */
    int   nc_nalloc;     /* Allocation size for the pointers above      */

    float *bbcz_simple;  /* bbcz for simple grid converted from super   */
    float *bb_simple;    /* bb for simple grid converted from super     */
    int   *flags_simple; /* flags for simple grid converted from super  */
    int   nc_nalloc_simple; /* Allocation size for the pointers above   */

    int  nsubc_tot;      /* Total number of subcell, used for printing  */
} nbnxn_grid_t;

#ifdef GMX_NBNXN_SIMD
#if GMX_NBNXN_SIMD_BITWIDTH == 128
#define GMX_MM128_HERE
#else
#if GMX_NBNXN_SIMD_BITWIDTH == 256
#define GMX_MM256_HERE
#else
#error "unsupported GMX_NBNXN_SIMD_BITWIDTH"
#endif
#endif
#include "gmx_simd_macros.h"

typedef struct nbnxn_x_ci_simd_4xn {
    /* The i-cluster coordinates for simple search */
    gmx_mm_pr ix_SSE0,iy_SSE0,iz_SSE0;
    gmx_mm_pr ix_SSE1,iy_SSE1,iz_SSE1;
    gmx_mm_pr ix_SSE2,iy_SSE2,iz_SSE2;
    gmx_mm_pr ix_SSE3,iy_SSE3,iz_SSE3;
} nbnxn_x_ci_simd_4xn_t;

typedef struct nbnxn_x_ci_simd_2xnn {
    /* The i-cluster coordinates for simple search */
    gmx_mm_pr ix_SSE0,iy_SSE0,iz_SSE0;
    gmx_mm_pr ix_SSE2,iy_SSE2,iz_SSE2;
} nbnxn_x_ci_simd_2xnn_t;

#endif

/* Working data for the actual i-supercell during pair search */
typedef struct nbnxn_list_work {
    gmx_cache_protect_t cp0; /* Protect cache between threads               */

    float *bb_ci;      /* The bounding boxes, pbc shifted, for each cluster */
    real  *x_ci;       /* The coordinates, pbc shifted, for each atom       */
#ifdef GMX_NBNXN_SIMD
    nbnxn_x_ci_simd_4xn_t *x_ci_simd_4xn;
    nbnxn_x_ci_simd_2xnn_t *x_ci_simd_2xnn;
#endif
    int  cj_ind;       /* The current cj_ind index for the current list     */
    int  cj4_init;     /* The first unitialized cj4 block                   */

    float *d2;         /* Bounding box distance work array                  */

    nbnxn_cj_t *cj;    /* The j-cell list                                   */
    int  cj_nalloc;    /* Allocation size of cj                             */

    int ncj_noq;       /* Nr. of cluster pairs without Coul for flop count  */
    int ncj_hlj;       /* Nr. of cluster pairs with 1/2 LJ for flop count   */

    gmx_cache_protect_t cp1; /* Protect cache between threads               */
} nbnxn_list_work_t;

/* Function type for setting the i-atom coordinate working data */
typedef void
gmx_icell_set_x_t(int ci,
                  real shx,real shy,real shz,
                  int na_c,
                  int stride,const real *x,
                  nbnxn_list_work_t *work);

static gmx_icell_set_x_t icell_set_x_simple;
#ifdef GMX_NBNXN_SIMD
static gmx_icell_set_x_t icell_set_x_simple_simd_4xn;
static gmx_icell_set_x_t icell_set_x_simple_simd_2xnn;
#endif
static gmx_icell_set_x_t icell_set_x_supersub;
#ifdef NBNXN_SEARCH_SSE
static gmx_icell_set_x_t icell_set_x_supersub_sse8;
#endif

#undef GMX_MM128_HERE
#undef GMX_MM256_HERE

/* Local cycle count struct for profiling */
typedef struct {
    int          count;
    gmx_cycles_t c;
    gmx_cycles_t start;
} nbnxn_cycle_t;

/* Local cycle count enum for profiling */
enum { enbsCCgrid, enbsCCsearch, enbsCCcombine, enbsCCreducef, enbsCCnr };

/* Thread-local work struct, contains part of nbnxn_grid_t */
typedef struct {
    gmx_cache_protect_t cp0;

    int *cxy_na;
    int cxy_na_nalloc;

    int  *sort_work;
    int  sort_work_nalloc;

    nbnxn_buffer_flags_t buffer_flags; /* Flags for force buffer access */

    int  ndistc;         /* Number of distance checks for flop counting */

    nbnxn_cycle_t cc[enbsCCnr];

    gmx_cache_protect_t cp1;
} nbnxn_search_work_t;

/* Main pair-search struct, contains the grid(s), not the pair-list(s) */
typedef struct nbnxn_search {
    int  ePBC;            /* PBC type enum                              */
    matrix box;           /* The periodic unit-cell                     */

    gmx_bool DomDec;      /* Are we doing domain decomposition?         */
    ivec dd_dim;          /* Are we doing DD in x,y,z?                  */
    gmx_domdec_zones_t *zones; /* The domain decomposition zones        */

    int  ngrid;           /* The number of grids, equal to #DD-zones    */
    nbnxn_grid_t *grid;   /* Array of grids, size ngrid                 */
    int  *cell;           /* Actual allocated cell array for all grids  */
    int  cell_nalloc;     /* Allocation size of cell                    */
    int  *a;              /* Atom index for grid, the inverse of cell   */
    int  a_nalloc;        /* Allocation size of a                       */

    int  natoms_local;    /* The local atoms run from 0 to natoms_local */
    int  natoms_nonlocal; /* The non-local atoms run from natoms_local
                           * to natoms_nonlocal */

    gmx_bool print_cycles;
    int      search_count;
    nbnxn_cycle_t cc[enbsCCnr];

    gmx_icell_set_x_t *icell_set_x; /* Function for setting i-coords    */

    int  nthread_max;     /* Maximum number of threads for pair-search  */
    nbnxn_search_work_t *work; /* Work array, size nthread_max          */
} nbnxn_search_t_t;


static void nbs_cycle_start(nbnxn_cycle_t *cc)
{
    cc->start = gmx_cycles_read();
}

static void nbs_cycle_stop(nbnxn_cycle_t *cc)
{
    cc->c += gmx_cycles_read() - cc->start;
    cc->count++;
}


#ifdef __cplusplus
}
#endif

#endif
