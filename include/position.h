/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
/*! \file
 * \brief API for handling positions.
 */
#ifndef POSITION_H
#define POSITION_H

#include "typedefs.h"

#include "indexutil.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Stores a set of positions together with their origins.
 */
typedef struct gmx_ana_pos_t
{
    /*! \brief
     * Number of positions.
     */
    int                 nr;
    /*! \brief
     * Array of positions.
     */
    rvec               *x;
    /*! \brief
     * Velocities (can be NULL).
     */
    rvec               *v;
    /*! \brief
     * Forces (can be NULL).
     */
    rvec               *f;
    /*! \brief
     * Mapping of the current positions to the original group.
     *
     * \see gmx_ana_indexmap_t
     */
    gmx_ana_indexmap_t  m;
    /*! \brief
     * Pointer to the current evaluation group.
     */
    gmx_ana_index_t    *g;
    /*! \brief
     * Number of elements allocated for \c x.
     */
    int                 nalloc_x;
} gmx_ana_pos_t;

/** Initializes an empty position structure. */
void
gmx_ana_pos_clear(gmx_ana_pos_t *pos);
/** Ensures that enough memory has been allocated to store positions. */
void
gmx_ana_pos_reserve(gmx_ana_pos_t *pos, int n, int isize);
/** Request memory allocation for velocities. */
void
gmx_ana_pos_reserve_velocities(gmx_ana_pos_t *pos);
/** Request memory allocation for forces. */
void
gmx_ana_pos_reserve_forces(gmx_ana_pos_t *pos);
/** Initializes a \c gmx_ana_pos_t to represent a constant position. */
void
gmx_ana_pos_init_const(gmx_ana_pos_t *pos, rvec x);
/** Frees the memory allocated for position storage. */
void
gmx_ana_pos_deinit(gmx_ana_pos_t *pos);
/** Frees the memory allocated for positions. */
void
gmx_ana_pos_free(gmx_ana_pos_t *pos);
/** Copies the evaluated positions to a preallocated data structure. */
void
gmx_ana_pos_copy(gmx_ana_pos_t *dest, gmx_ana_pos_t *src, gmx_bool bFirst);

/** Sets the number of positions in a position structure. */
void
gmx_ana_pos_set_nr(gmx_ana_pos_t *pos, int n);
/** Sets the evaluation group of a position data structure. */
void
gmx_ana_pos_set_evalgrp(gmx_ana_pos_t *pos, gmx_ana_index_t *g);
/** Empties a position data structure with full initialization. */
void
gmx_ana_pos_empty_init(gmx_ana_pos_t *pos);
/** Empties a position data structure. */
void
gmx_ana_pos_empty(gmx_ana_pos_t *pos);
/** Appends a position to a preallocated data structure with full
 * initialization. */
void
gmx_ana_pos_append_init(gmx_ana_pos_t *dest, gmx_ana_index_t *g,
                        gmx_ana_pos_t *src, int i);
/** Appends a position to a preallocated data structure. */
void
gmx_ana_pos_append(gmx_ana_pos_t *dest, gmx_ana_index_t *g,
                   gmx_ana_pos_t *src, int i, int refid);
/** Updates position data structure state after appends. */
void
gmx_ana_pos_append_finish(gmx_ana_pos_t *pos);

#ifdef __cplusplus
}
#endif

#endif
