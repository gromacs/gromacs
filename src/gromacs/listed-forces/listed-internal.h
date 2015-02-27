/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/bitmask.h"

/*! \brief struct with all output for bonded forces, used per thread */
typedef struct
{
    rvec             *f;            /**< Force array */
    int               f_nalloc;     /**< Allocation size of f */
    gmx_bitmask_t     red_mask;     /**< Mask for marking which parts of f are filled */
    rvec             *fshift;       /**< Shift force array, size SHIFTS */
    real              ener[F_NRE];  /**< Energy array */
    gmx_grppairener_t grpp;         /**< Group pair energy data for pairs */
    real              dvdl[efptNR]; /**< Free-energy dV/dl output */
}
f_thread_t;

/*! \brief struct contain all data for bonded force threading */
struct bonded_threading_t
{
    /* Thread local force and energy data */
    int         nthreads;   /**< Number of threads to be used for bondeds */
    int         red_ashift; /**< Size of force reduction blocks in bits */
    int         red_nblock; /**< The number of force blocks to reduce */
    f_thread_t *f_t;        /**< Force/enegry data per thread, size nthreads */

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
