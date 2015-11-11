/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

/* Struct for passing around the nbnxn grid corners and atom density */
struct nbnxn_grid_corners_density_t {
    rvec corner0;
    rvec corner1;
    real atom_density;
};

/* Struct for passing around all grid dimensions */
struct nbnxn_grid_dims_t {
    int  ncolumn_x;    /**< Grid size along x */
    int  ncolumn_y;    /**< Grid size along y */
    rvec corner0;      /**< Lower grid corner */
    rvec corner1;      /**< Upper grid corner */
    real size_x;       /**< Grid cell sizes along x */
    real size_y;       /**< Grid cell sizes along y */
    int  natoms;       /**< The number of atoms to communicate */
    ivec work;         /**< Storage space used to communicate extra data during the DD halo communication setup */
};

/* Allocate and initialize the pair search grids for all zones in nbs */
void nbnxn_grids_init(nbnxn_search_t nbs, int nzone);

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
                       int grid_ind,
                       const nbnxn_grid_corners_density_t *corners_density,
                       int a0, int a1,
                       const int *atinfo,
                       rvec *x,
                       int nmoved, int *move,
                       int nb_kernel_type,
                       real pairlist_cutoff,
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

/* Add simple grid type information to the local super/sub grid */
void nbnxn_grid_add_simple(nbnxn_search_t    nbs,
                           nbnxn_atomdata_t *nbat);

/* Return the dimensions of the local grid in dims */
void nbnxn_get_local_grid_dimensions(nbnxn_search_t     nbs,
                                     nbnxn_grid_dims_t *dims);

/* Struct for returning all information about a single grid column */
struct nbnxn_grid_column_t {
    nbnxn_bb_t        column_bb;   /* Bounding box of the column */
    int               bb_start;    /* Index of the first bounding box / cluster */
    int               nbb;         /* Number of bounding boxes in bb */
    const nbnxn_bb_t *bb;          /* Bounding boxes in this column are bb[0] - bb[nbb-1], is NULL with GPU grids */
    const float      *bb_z;        /* Bounding boxes z-coordinates in this column are bbz[0] - bbz[nbb], is NULL with CPU grids */
    int               atom_start;  /* First nbnxn atom index in this column */
    int               bb_natoms;   /* Number of atoms per bounding box */
    int               natoms;      /* Number of real atoms (not fillers) in this column */
};

/* Return the dimensions and contents of grid column with index cx, cy
 * of the local pair search grid.
 */
void nbnxn_get_local_grid_column(nbnxn_search_t nbs, int cx, int cy,
                                 nbnxn_grid_column_t *col);

/* Add a grid for DD zone dd_zone to nbs and set the grid dimensions and
 * atom count for each column.
 * Note that this does not set the atom indices (not required for non-local)
 * and the bounding boxes (need to be calculated later).
 */
void nbnxn_set_grid_parameters(nbnxn_search_t           nbs,
                               int                      dd_zone,
                               const nbnxn_grid_dims_t *dims,
                               int                      atom_start,
                               const int               *column_natoms);

/* Return the order indices *a of the atoms on the ns grid, size n */
void nbnxn_get_atomorder(const nbnxn_search_t nbs, const int **a, int *n);

/* Renumber the atom indices on the grid to consecutive order */
void nbnxn_set_atomorder(nbnxn_search_t nbs);

#endif
