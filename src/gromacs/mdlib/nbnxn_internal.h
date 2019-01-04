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

#ifndef _nbnxn_internal_h
#define _nbnxn_internal_h

#include <memory>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_zones_t;


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
template<int packSize> static inline int atom_to_x_index(int a)
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


template <class T>
using AlignedVector = std::vector < T, gmx::AlignedAllocator < T>>;


/* Bounding box for a nbnxn atom cluster */
typedef struct {
    float lower[NNBSBB_C];
    float upper[NNBSBB_C];
} nbnxn_bb_t;


/* A pair-search grid struct for one domain decomposition zone
 *
 * Note that when atom groups, instead of individual atoms, are assigned
 * to grid cells, individual atoms can be geometrically outside the cell
 * and grid that they have been assigned to (as determined by the center
 * or geometry of the atom group they belong to).
 */
struct nbnxn_grid_t
{
    rvec     c0;                   /* The lower corner of the (local) grid        */
    rvec     c1;                   /* The upper corner of the (local) grid        */
    rvec     size;                 /* c1 - c0                                     */
    real     atom_density;         /* The atom number density for the local grid  */
    real     maxAtomGroupRadius;   /* The maximum distance an atom can be outside
                                    * of a cell and outside of the grid
                                    */

    gmx_bool bSimple;              /* Is this grid simple or super/sub            */
    int      na_c;                 /* Number of atoms per cluster                 */
    int      na_cj;                /* Number of atoms for list j-clusters         */
    int      na_sc;                /* Number of atoms per super-cluster           */
    int      na_c_2log;            /* 2log of na_c                                */

    int      numCells[DIM - 1];    /* Number of cells along x/y                   */
    int      nc;                   /* Total number of cells                       */

    real     cellSize[DIM - 1];    /* size of a cell                              */
    real     invCellSize[DIM - 1]; /* 1/cellSize                                  */

    int      cell0;                /* Index in nbs->cell corresponding to cell 0  */

    /* Grid data */
    std::vector<int> cxy_na;        /* The number of atoms for each column in x,y  */
    std::vector<int> cxy_ind;       /* Grid (super)cell index, offset from cell0   */

    std::vector<int> nsubc;         /* The number of sub cells for each super cell */

    /* Bounding boxes */
    std::vector<float>                                    bbcz;                /* Bounding boxes in z for the cells */
    std::vector < nbnxn_bb_t, gmx::AlignedAllocator < nbnxn_bb_t>> bb;         /* 3D bounding boxes for the sub cells */
    std::vector < nbnxn_bb_t, gmx::AlignedAllocator < nbnxn_bb_t>> bbjStorage; /* 3D j-bounding boxes for the case where
                                                                                * the i- and j-cluster sizes are different */
    gmx::ArrayRef<nbnxn_bb_t>                              bbj;                /* 3D j-bounding boxes */
    std::vector < float, gmx::AlignedAllocator < float>>            pbb;       /* 3D b. boxes in xxxx format per super cell   */

    /* Bit-flag information */
    std::vector<int>          flags;     /* Flags for properties of clusters in each cell */
    std::vector<unsigned int> fep;       /* FEP signal bits for atoms in each cluster */

    /* Statistics */
    int                       nsubc_tot; /* Total number of subcell, used for printing  */
};

/* Working data for the actual i-supercell during pair search */
struct NbnxnPairlistCpuWork
{
    // Struct for storing coordinats and bounding box for an i-entry during search
    struct IClusterData
    {
        IClusterData() :
            bb(1),
            x(c_nbnxnCpuIClusterSize*DIM),
            xSimd(c_nbnxnCpuIClusterSize*DIM*GMX_REAL_MAX_SIMD_WIDTH)
        {
        }

        // The bounding boxes, pbc shifted, for each cluster
        AlignedVector<nbnxn_bb_t> bb;
        // The coordinates, pbc shifted, for each atom
        std::vector<real>         x;
        // Aligned list for storing 4*DIM*GMX_SIMD_REAL_WIDTH reals
        AlignedVector<real>       xSimd;
    };

    // Protect data from cache pollution between threads
    gmx_cache_protect_t       cp0;

    // Work data for generating an IEntry in the pairlist
    IClusterData              iClusterData;
    // The current cj_ind index for the current list
    int                       cj_ind;
    // Temporary j-cluster list, used for sorting on exclusions
    std::vector<nbnxn_cj_t>   cj;

    // Nr. of cluster pairs without Coulomb for flop counting
    int                       ncj_noq;
    // Nr. of cluster pairs with 1/2 LJ for flop count
    int                       ncj_hlj;

    // Protect data from cache pollution between threads
    gmx_cache_protect_t       cp1;
};

/* Working data for the actual i-supercell during pair search */
struct NbnxnPairlistGpuWork
{
    struct ISuperClusterData
    {
        ISuperClusterData() :
            bb(c_gpuNumClusterPerCell),
#if NBNXN_SEARCH_BB_SIMD4
            bbPacked(c_gpuNumClusterPerCell/STRIDE_PBB*NNBSBB_XXXX),
#endif
            x(c_gpuNumClusterPerCell*c_nbnxnGpuClusterSize*DIM),
            xSimd(c_gpuNumClusterPerCell*c_nbnxnGpuClusterSize*DIM)
        {
        }

        // The bounding boxes, pbc shifted, for each cluster
        AlignedVector<nbnxn_bb_t> bb;
        // As bb, but in packed xxxx format
        AlignedVector<float>      bbPacked;
        // The coordinates, pbc shifted, for each atom
        AlignedVector<real>       x;
        // Aligned coordinate list used for 4*DIM*GMX_SIMD_REAL_WIDTH floats
        AlignedVector<real>       xSimd;
    };

    NbnxnPairlistGpuWork() :
        distanceBuffer(c_gpuNumClusterPerCell),
        sci_sort({}, {gmx::PinningPolicy::PinnedIfSupported})
    {
    }

    // Protect data from cache pollution between threads
    gmx_cache_protect_t       cp0;

    // Work data for generating an i-entry in the pairlist
    ISuperClusterData         iSuperClusterData;
    // The current j-cluster index for the current list
    int                       cj_ind;
    // Bounding box distance work array
    AlignedVector<float>      distanceBuffer;

    // Buffer for sorting list entries
    std::vector<int>          sortBuffer;

    // Second sci array, for sorting
    gmx::HostVector<nbnxn_sci_t> sci_sort;

    // Protect data from cache pollution between threads
    gmx_cache_protect_t       cp1;
};

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
struct nbnxn_search_work_t
{
    nbnxn_search_work_t();

    ~nbnxn_search_work_t();

    gmx_cache_protect_t       cp0;          /* Buffer to avoid cache polution */

    std::vector<int>          cxy_na;       /* Grid column atom counts temporary buffer */

    std::vector<int>          sortBuffer;   /* Temporary buffer for sorting atoms within a grid column */

    nbnxn_buffer_flags_t      buffer_flags; /* Flags for force buffer access */

    int                       ndistc;       /* Number of distance checks for flop counting */


    std::unique_ptr<t_nblist> nbl_fep;      /* Temporary FEP list for load balancing */

    nbnxn_cycle_t             cc[enbsCCnr]; /* Counters for thread-local cycles */

    gmx_cache_protect_t       cp1;          /* Buffer to avoid cache polution */
};

/* Main pair-search struct, contains the grid(s), not the pair-list(s) */
struct nbnxn_search
{
    /* \brief Constructor
     *
     * \param[in] n_dd_cells   The number of domain decomposition cells per dimension, without DD nullptr should be passed
     * \param[in] zones        The domain decomposition zone setup, without DD nullptr should be passed
     * \param[in] bFEP         Tells whether non-bonded interactions are perturbed
     * \param[in] nthread_max  The maximum number of threads used in the search
     */

    nbnxn_search(const ivec               *n_dd_cells,
                 const gmx_domdec_zones_t *zones,
                 gmx_bool                  bFEP,
                 int                       nthread_max);

    gmx_bool                   bFEP;            /* Do we have perturbed atoms? */
    int                        ePBC;            /* PBC type enum                              */
    matrix                     box;             /* The periodic unit-cell                     */

    gmx_bool                   DomDec;          /* Are we doing domain decomposition?         */
    ivec                       dd_dim;          /* Are we doing DD in x,y,z?                  */
    const gmx_domdec_zones_t  *zones;           /* The domain decomposition zones        */

    std::vector<nbnxn_grid_t>  grid;            /* Array of grids, size ngrid                 */
    std::vector<int>           cell;            /* Actual allocated cell array for all grids  */
    std::vector<int>           a;               /* Atom index for grid, the inverse of cell   */

    int                        natoms_local;    /* The local atoms run from 0 to natoms_local */
    int                        natoms_nonlocal; /* The non-local atoms run from natoms_local
                                                 * to natoms_nonlocal */

    gmx_bool             print_cycles;
    int                  search_count;
    nbnxn_cycle_t        cc[enbsCCnr];

    /* Thread-local work data */
    mutable std::vector<nbnxn_search_work_t> work; /* Work array, one entry for each thread */
};


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
