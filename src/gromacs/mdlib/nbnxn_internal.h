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

#ifndef _nbnxn_internal_h
#define _nbnxn_internal_h

#include "gromacs/domdec/domdec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/real.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

struct gmx_domdec_zones_t;


/* The number of clusters in a pair-search cell, used for GPU */
static const int c_gpuNumClusterPerCellZ = 2;
static const int c_gpuNumClusterPerCellY = 2;
static const int c_gpuNumClusterPerCellX = 2;
static const int c_gpuNumClusterPerCell  = c_gpuNumClusterPerCellZ*c_gpuNumClusterPerCellY*c_gpuNumClusterPerCellX;


/* Strides for x/f with xyz and xyzq coordinate (and charge) storage */
#define STRIDE_XYZ         3
#define STRIDE_XYZQ        4
/* Size of packs of x, y or z with SIMD packed coords/forces */
static const int c_packX4 = 4;
static const int c_packX8 = 8;
/* Strides for a pack of 4 and 8 coordinates/forces */
#define STRIDE_P4         (DIM*c_packX4)
#define STRIDE_P8         (DIM*c_packX8)

/* Returns the index in a coordinate array corresponding to atom a */
template<int packSize> static gmx_inline int atom_to_x_index(int a)
{
    return DIM*(a & ~(packSize - 1)) + (a & (packSize - 1));
}


#if GMX_SIMD
/* Memory alignment in bytes as required by SIMD aligned loads/stores */
#define NBNXN_MEM_ALIGN  (GMX_SIMD_REAL_WIDTH*sizeof(real))
#else
/* No alignment required, but set it so we can call the same routines */
#define NBNXN_MEM_ALIGN  32
#endif


/* Bounding box calculations are (currently) always in single precision, so
 * we only need to check for single precision support here.
 * This uses less (cache-)memory and SIMD is faster, at least on x86.
 */
#if GMX_SIMD4_HAVE_FLOAT
#    define NBNXN_SEARCH_BB_SIMD4      1
/* Memory alignment in bytes as required by SIMD aligned loads/stores */
#    define NBNXN_SEARCH_BB_MEM_ALIGN  (GMX_SIMD4_WIDTH*sizeof(float))
#else
#    define NBNXN_SEARCH_BB_SIMD4      0
/* No alignment required, but set it so we can call the same routines */
#    define NBNXN_SEARCH_BB_MEM_ALIGN  32
#endif


/* Pair search box lower and upper corner in x,y,z.
 * Store this in 4 iso 3 reals, which is useful with 4-wide SIMD.
 * To avoid complicating the code we also use 4 without 4-wide SIMD.
 */
#define NNBSBB_C         4
/* Pair search box lower and upper bound in z only. */
#define NNBSBB_D         2
/* Pair search box lower and upper corner x,y,z indices, entry 3 is unused */
#define BB_X  0
#define BB_Y  1
#define BB_Z  2


#if NBNXN_SEARCH_BB_SIMD4
/* Always use 4-wide SIMD for bounding box calculations */

#    if !GMX_DOUBLE
/* Single precision BBs + coordinates, we can also load coordinates with SIMD */
#        define NBNXN_SEARCH_SIMD4_FLOAT_X_BB  1
#    else
#        define NBNXN_SEARCH_SIMD4_FLOAT_X_BB  0
#    endif

/* The packed bounding box coordinate stride is always set to 4.
 * With AVX we could use 8, but that turns out not to be faster.
 */
#    define STRIDE_PBB       GMX_SIMD4_WIDTH
#    define STRIDE_PBB_2LOG  2

/* Store bounding boxes corners as quadruplets: xxxxyyyyzzzz */
#    define NBNXN_BBXXXX  1
/* Size of bounding box corners quadruplet */
#    define NNBSBB_XXXX  (NNBSBB_D*DIM*STRIDE_PBB)

#else  /* NBNXN_SEARCH_BB_SIMD4 */

#    define NBNXN_SEARCH_SIMD4_FLOAT_X_BB  0
#    define NBNXN_BBXXXX                   0

#endif /* NBNXN_SEARCH_BB_SIMD4 */


/* Bounding box for a nbnxn atom cluster */
typedef struct {
    float lower[NNBSBB_C];
    float upper[NNBSBB_C];
} nbnxn_bb_t;


/* A pair-search grid struct for one domain decomposition zone */
typedef struct {
    rvec          c0;               /* The lower corner of the (local) grid        */
    rvec          c1;               /* The upper corner of the (local) grid        */
    rvec          size;             /* c1 - c0                                     */
    real          atom_density;     /* The atom number density for the local grid  */

    gmx_bool      bSimple;          /* Is this grid simple or super/sub            */
    int           na_c;             /* Number of atoms per cluster                 */
    int           na_cj;            /* Number of atoms for list j-clusters         */
    int           na_sc;            /* Number of atoms per super-cluster           */
    int           na_c_2log;        /* 2log of na_c                                */

    int           ncx;              /* Number of (super-)cells along x             */
    int           ncy;              /* Number of (super-)cells along y             */
    int           nc;               /* Total number of (super-)cells               */

    real          sx;               /* x-size of a (super-)cell                    */
    real          sy;               /* y-size of a (super-)cell                    */
    real          inv_sx;           /* 1/sx                                        */
    real          inv_sy;           /* 1/sy                                        */

    int           cell0;            /* Index in nbs->cell corresponding to cell 0  */

    int          *cxy_na;           /* The number of atoms for each column in x,y  */
    int          *cxy_ind;          /* Grid (super)cell index, offset from cell0   */
    int           cxy_nalloc;       /* Allocation size for cxy_na and cxy_ind      */

    int          *nsubc;            /* The number of sub cells for each super cell */
    float        *bbcz;             /* Bounding boxes in z for the super cells     */
    nbnxn_bb_t   *bb;               /* 3D bounding boxes for the sub cells         */
    nbnxn_bb_t   *bbj;              /* 3D j-bounding boxes for the case where      *
                                     * the i- and j-cluster sizes are different    */
    float        *pbb;              /* 3D b. boxes in xxxx format per super cell   */
    int          *flags;            /* Flag for the super cells                    */
    unsigned int *fep;              /* FEP signal bits for sub cells               */
    int           nc_nalloc;        /* Allocation size for the pointers above      */

    int           nsubc_tot;        /* Total number of subcell, used for printing  */
} nbnxn_grid_t;

/* Working data for the actual i-supercell during pair search */
typedef struct nbnxn_list_work {
    gmx_cache_protect_t     cp0;             /* Protect cache between threads               */

    nbnxn_bb_t             *bb_ci;           /* The bounding boxes, pbc shifted, for each cluster */
    float                  *pbb_ci;          /* As bb_ci, but in xxxx packed format               */
    real                   *x_ci;            /* The coordinates, pbc shifted, for each atom       */
    real                   *x_ci_simd;       /* aligned pointer to 4*DIM*GMX_SIMD_REAL_WIDTH floats */
    int                     cj_ind;          /* The current cj_ind index for the current list     */
    int                     cj4_init;        /* The first unitialized cj4 block                   */

    float                  *d2;              /* Bounding box distance work array                  */

    nbnxn_cj_t             *cj;              /* The j-cell list                                   */
    int                     cj_nalloc;       /* Allocation size of cj                             */

    int                     ncj_noq;         /* Nr. of cluster pairs without Coul for flop count  */
    int                     ncj_hlj;         /* Nr. of cluster pairs with 1/2 LJ for flop count   */

    int                    *sort;            /* Sort index                    */
    int                     sort_nalloc;     /* Allocation size of sort       */

    nbnxn_sci_t            *sci_sort;        /* Second sci array, for sorting */
    int                     sci_sort_nalloc; /* Allocation size of sci_sort   */

    gmx_cache_protect_t     cp1;             /* Protect cache between threads               */
} nbnxn_list_work_t;

/* Function type for setting the i-atom coordinate working data */
typedef void
    gmx_icell_set_x_t (int ci,
                       real shx, real shy, real shz,
                       int stride, const real *x,
                       nbnxn_list_work_t *work);

/* Local cycle count struct for profiling */
typedef struct {
    int          count;
    gmx_cycles_t c;
    gmx_cycles_t start;
} nbnxn_cycle_t;

/* Local cycle count enum for profiling */
enum {
    enbsCCgrid, enbsCCsearch, enbsCCcombine, enbsCCreducef, enbsCCnr
};

/* Thread-local work struct, contains part of nbnxn_grid_t */
typedef struct {
    gmx_cache_protect_t  cp0;

    int                 *cxy_na;
    int                  cxy_na_nalloc;

    int                 *sort_work;
    int                  sort_work_nalloc;

    nbnxn_buffer_flags_t buffer_flags; /* Flags for force buffer access */

    int                  ndistc;       /* Number of distance checks for flop counting */

    t_nblist            *nbl_fep;      /* Temporary FEP list for load balancing */

    nbnxn_cycle_t        cc[enbsCCnr];

    gmx_cache_protect_t  cp1;
} nbnxn_search_work_t;

/* Main pair-search struct, contains the grid(s), not the pair-list(s) */
typedef struct nbnxn_search {
    gmx_bool                   bFEP;            /* Do we have perturbed atoms? */
    int                        ePBC;            /* PBC type enum                              */
    matrix                     box;             /* The periodic unit-cell                     */

    gmx_bool                   DomDec;          /* Are we doing domain decomposition?         */
    ivec                       dd_dim;          /* Are we doing DD in x,y,z?                  */
    struct gmx_domdec_zones_t *zones;           /* The domain decomposition zones        */

    int                        ngrid;           /* The number of grids, equal to #DD-zones    */
    nbnxn_grid_t              *grid;            /* Array of grids, size ngrid                 */
    int                       *cell;            /* Actual allocated cell array for all grids  */
    int                        cell_nalloc;     /* Allocation size of cell                    */
    int                       *a;               /* Atom index for grid, the inverse of cell   */
    int                        a_nalloc;        /* Allocation size of a                       */

    int                        natoms_local;    /* The local atoms run from 0 to natoms_local */
    int                        natoms_nonlocal; /* The non-local atoms run from natoms_local
                                                 * to natoms_nonlocal */

    gmx_bool             print_cycles;
    int                  search_count;
    nbnxn_cycle_t        cc[enbsCCnr];

    gmx_icell_set_x_t   *icell_set_x; /* Function for setting i-coords    */

    int                  nthread_max; /* Maximum number of threads for pair-search  */
    nbnxn_search_work_t *work;        /* Work array, size nthread_max          */
} nbnxn_search_t_t;


static inline void nbs_cycle_start(nbnxn_cycle_t *cc)
{
    cc->start = gmx_cycles_read();
}

static inline void nbs_cycle_stop(nbnxn_cycle_t *cc)
{
    cc->c += gmx_cycles_read() - cc->start;
    cc->count++;
}


#endif
