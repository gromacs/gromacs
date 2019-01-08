/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

/*! \internal \file
 *
 * \brief Declares the grid and bounding box objects
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_nbnxn
 */

#ifndef GMX_NBNXN_GRID_H
#define GMX_NBNXN_GRID_H

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"


struct gmx_domdec_zones_t;


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


/* Bounding box for a nbnxn atom cluster */
typedef struct {
    float lower[NNBSBB_C];
    float upper[NNBSBB_C];
} nbnxn_bb_t;


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

#endif
