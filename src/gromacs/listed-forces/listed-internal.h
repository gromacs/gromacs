/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * \brief This file contains declarations for functions needed
 * internally by the module.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_listed-forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_INTERNAL_H
#define GMX_LISTED_FORCES_LISTED_INTERNAL_H

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/bitmask.h"

/* We reduce the force array in blocks of 32 atoms. This is large enough
 * to not cause overhead and 32*sizeof(rvec) is a multiple of the cache-line
 * size on all systems.
 */
static const int reduction_block_size = 32; /**< Force buffer block size in atoms*/
static const int reduction_block_bits =  5; /**< log2(reduction_block_size) */

/*! \internal \brief struct with output for bonded forces, used per thread */
typedef struct
{
    rvec4            *f;            /**< Force array */
    int               f_nalloc;     /**< Allocation size of f */
    gmx_bitmask_t    *mask;         /**< Mask for marking which parts of f are filled, working array for constructing mask in bonded_threading_t */
    int               nblock_used;  /**< Number of blocks touched by our thread */
    int              *block_index;  /**< Index to touched blocks, size nblock_used */
    int               block_nalloc; /**< Allocation size of f (*reduction_block_size), mask_index, mask */

    rvec             *fshift;       /**< Shift force array, size SHIFTS */
    real              ener[F_NRE];  /**< Energy array */
    gmx_grppairener_t grpp;         /**< Group pair energy data for pairs */
    real              dvdl[efptNR]; /**< Free-energy dV/dl output */
}
f_thread_t;

/*! \internal \brief struct contain all data for bonded force threading */
struct bonded_threading_t
{
    /* Thread local force and energy data */
    int            nthreads;     /**< Number of threads to be used for bondeds */
    f_thread_t    *f_t;          /**< Force/enegry data per thread, size nthreads */
    int            nblock_used;  /**< The number of force blocks to reduce */
    int           *block_index;  /**< Index of size nblock_used into mask */
    gmx_bitmask_t *mask;         /**< Mask array, one element corresponds to a block of reduction_block_size atoms of the force array, bit corresponding to thread indices set if a thread writes to that block */
    int            block_nalloc; /**< Allocation size of block_index and mask */

    /* There are two different ways to distribute the bonded force calculation
     * over the threads. We dedice which to use based on the number of threads.
     */
    int bonded_max_nthread_uniform; /**< Maximum thread count for uniform distribution of bondeds over threads */
};


/*! \brief Returns the global topology atom number belonging to local
 * atom index i.
 *
 * This function is intended for writing ascii output and returns atom
 * numbers starting at 1.  When global_atom_index=NULL returns i+1.
 */
int
glatnr(int *global_atom_index, int i);

#endif
