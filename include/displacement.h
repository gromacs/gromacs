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
 * \brief API for on-line calculation of displacements.
 *
 * The API is documented in more detail on a separate page:
 * \ref displacements
 *
 * The functions within this file can be used and developed independently of
 * the other parts of the library.
 * Other parts of the library do not reference these functions.
 */
#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Data structure for displacement calculation. */
typedef struct gmx_ana_displ_t gmx_ana_displ_t;

struct gmx_ana_pos_t;

/** Allocates memory for displacement calculation. */
int
gmx_ana_displ_create(gmx_ana_displ_t **d, int nmax, real tmax);
/** Initializes displacement calculation for a frame. */
int
gmx_ana_displ_start_frame(gmx_ana_displ_t *d, real t);
/** Returns the number of steps corresponding to a given time interval. */
int
gmx_ana_displ_time_to_steps(gmx_ana_displ_t *d, real time, int *nsteps);
/** Stores the position of a particle for displacement calculation. */
int
gmx_ana_displ_store(gmx_ana_displ_t *d, atom_id id, rvec x, gmx_bool bPres);
/** Convenience function for storing an array of particle positions for displacement calculation. */
int
gmx_ana_displ_store_array(gmx_ana_displ_t *d, int n, atom_id id[], rvec x[]);
/** Stores an array of particle positions for displacement calculation, including unselected particles. */
int
gmx_ana_displ_store_all(gmx_ana_displ_t *d, atom_id id[], rvec x[]);
/** Convenience function for storing a set of positions from \c gmx_ana_pos_t. */
int
gmx_ana_displ_store_pos(gmx_ana_displ_t *d, struct gmx_ana_pos_t *p);
/** Calculates the displacement vector for a particle. */
int
gmx_ana_displ_vector(gmx_ana_displ_t *d, int step, t_pbc *pbc,
                     atom_id id, rvec x, rvec xout, gmx_bool *pout);
/** Calculates the displacement vectors for a list of particles. */
int
gmx_ana_displ_vectors(gmx_ana_displ_t *d, int step, t_pbc *pbc,
                      int n, atom_id id[], rvec x[],
                      rvec xout[], gmx_bool *pout);
/** Calculates the displacement vectors for all particles, including unselected. */
int
gmx_ana_displ_vectors_all(gmx_ana_displ_t *d, int step, t_pbc *pbc,
                          rvec x[], rvec xout[], gmx_bool *pout);
/** Frees the memory allocated for displacement calculation. */
void
gmx_ana_displ_free(gmx_ana_displ_t *d);

#ifdef __cplusplus
}
#endif

#endif
