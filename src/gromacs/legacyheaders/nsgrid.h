/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Used when estimating the interaction density.
 *
 * GRID_STDDEV_FAC * stddev estimates the interaction density. The
 * value sqrt(3) == 1.73205080757 gives a uniform load for a
 * rectangular 3D block of charge groups. For a sphere, it is not a
 * bad approximation for 4x1x1 up to 4x2x2.
 *
 * \todo It would be nicer to use sqrt(3) here, when all code that
 * includes this file is in C++, which will let us cope with the
 * std::sqrt<T> on Windows. */
static const real GRID_STDDEV_FAC = 1.73205080757;

/*! \brief The extent of the neighborsearch grid is a bit larger than sqrt(3)
 * to account for less dense regions at the edges of the system.
 */
static const real NSGRID_STDDEV_FAC = 2.0;

#define NSGRID_SIGNAL_MOVED_FAC  4
/* A cell index of NSGRID_SIGNAL_MOVED_FAC*ncells signals
 * that a charge group moved to another DD domain.
 */

t_grid *init_grid(FILE *fplog, t_forcerec *fr);

void done_grid(t_grid *grid);

void get_nsgrid_boundaries(int nboundeddim, matrix box,
                           gmx_domdec_t *dd,
                           gmx_ddbox_t *ddbox,
                           rvec *gr0, rvec *gr1,
                           int ncg, rvec *cgcm,
                           rvec grid_x0, rvec grid_x1,
                           real *grid_density);
/* Return the ns grid boundaries grid_x0 and grid_x1
 * and the estimate for the grid density.
 * For non-bounded dimensions the boundaries are determined
 * from the average and std.dev. of cgcm.
 * The are determined from box, unless gr0!=NULL or gr1!=NULL,
 * then they are taken from gr0 or gr1.
 * With dd and unbounded dimensions, the proper grid borders for cells
 * on the edges are determined from cgcm.
 */

void grid_first(FILE *log, t_grid *grid,
                gmx_domdec_t *dd, const gmx_ddbox_t *ddbox, matrix box, rvec izones_x0, rvec izones_x1,
                real rlong, real grid_density);

void fill_grid(gmx_domdec_zones_t *dd_zones,
               t_grid *grid, int ncg_tot,
               int cg0, int cg1, rvec cg_cm[]);
/* Allocates space on the grid for ncg_tot cg's.
 * Fills the grid with cg's from cg0 to cg1.
 * When cg0 is -1, contiues filling from grid->nr to cg1.
 */

void calc_elemnr(t_grid *grid, int cg0, int cg1, int ncg);

void calc_ptrs(t_grid *grid);

void grid_last(t_grid *grid, int cg0, int cg1, int ncg);

int xyz2ci_(int nry, int nrz, int x, int y, int z);
#define xyz2ci(nry, nrz, x, y, z) ((nry)*(nrz)*(x)+(nrz)*(y)+(z))
/* Return the cell index */

void ci2xyz(t_grid *grid, int i, int *x, int *y, int *z);

void check_grid(t_grid *grid);

void print_grid(FILE *log, t_grid *grid);

#ifdef __cplusplus
}
#endif
