/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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

#ifndef _nbnxn_grid_h
#define _nbnxn_grid_h

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_zones_t;

/* Allocate and initialize ngrid pair search grids in nbs */
void nbnxn_grids_init(nbnxn_search_t nbs, int ngrid);

/* Put the atoms on the pair search grid.
 * Only atoms a0 to a1 in x are put on the grid.
 * The atom_density is used to determine the grid size.
 * When atom_density<=0, the density is determined from a1-a0 and the corners.
 * With domain decomposition part of the n particles might have migrated,
 * but have not been removed yet. This count is given by nmoved.
 * When move[i] < 0 particle i has migrated and will not be put on the grid.
 * Without domain decomposition move will be NULL.
 */
void nbnxn_put_on_grid(nbnxn_search_t nbs,
                       int ePBC, matrix box,
                       int dd_zone,
                       rvec corner0, rvec corner1,
                       int a0, int a1,
                       real atom_density,
                       const int *atinfo,
                       rvec *x,
                       int nmoved, int *move,
                       int nb_kernel_type,
                       nbnxn_atomdata_t *nbat);

/* As nbnxn_put_on_grid, but for the non-local atoms
 * with domain decomposition. Should be called after calling
 * nbnxn_search_put_on_grid for the local atoms / home zone.
 */
void nbnxn_put_on_grid_nonlocal(nbnxn_search_t                   nbs,
                                const struct gmx_domdec_zones_t *zones,
                                const int                       *atinfo,
                                rvec                            *x,
                                int                              nb_kernel_type,
                                nbnxn_atomdata_t                *nbat);

/* Return the number of x and y cells in the local grid */
void nbnxn_get_ncells(nbnxn_search_t nbs, int *ncx, int *ncy);

/* Return the order indices *a of the atoms on the ns grid, size n */
void nbnxn_get_atomorder(const nbnxn_search_t nbs, const int **a, int *n);

/* Renumber the atom indices on the grid to consecutive order */
void nbnxn_set_atomorder(nbnxn_search_t nbs);

#endif
