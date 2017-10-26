/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
/* TODO find out what this file should be called */
#ifndef GMX_EWALD_PME_GRID_H
#define GMX_EWALD_PME_GRID_H

#include "config.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_pme_t;

/*! \brief
 * We allow coordinates to be out the unit-cell by up to 2 box lengths,
 * which might be needed along dimension x for a very skewed unit-cell.
 */
constexpr int c_pmeMaxUnitcellShift = 2;

/*! \brief
 * This affects the size of the lookup table of the modulo operation result,
 * when working with PME local grid indices of the particles.
 */
constexpr int c_pmeNeighborUnitcellCount = 2*c_pmeMaxUnitcellShift + 1;

struct pmegrid_t;
struct pmegrids_t;

#if GMX_MPI
void
gmx_sum_qgrid_dd(gmx_pme_t *pme, real *grid, int direction);
#endif

int
copy_pmegrid_to_fftgrid(const gmx_pme_t *pme, real *pmegrid, real *fftgrid, int grid_index);

int
copy_fftgrid_to_pmegrid(gmx_pme_t *pme, const real *fftgrid, real *pmegrid, int grid_index,
                        int nthread, int thread);

void
wrap_periodic_pmegrid(const gmx_pme_t *pme, real *pmegrid);

void
unwrap_periodic_pmegrid(gmx_pme_t *pme, real *pmegrid);

void
pmegrid_init(pmegrid_t *grid,
             int cx, int cy, int cz,
             int x0, int y0, int z0,
             int x1, int y1, int z1,
             gmx_bool set_alignment,
             int pme_order,
             real *ptr);

void
pmegrids_init(pmegrids_t *grids,
              int nx, int ny, int nz, int nz_base,
              int pme_order,
              gmx_bool bUseThreads,
              int nthread,
              int overlap_x,
              int overlap_y);

void
pmegrids_destroy(pmegrids_t *grids);

void
make_gridindex_to_localindex(int n, int local_start, int local_range,
                             int **global_to_local,
                             real **fraction_shift);

void
set_grid_alignment(int *pmegrid_nz, int pme_order);

void
reuse_pmegrids(const pmegrids_t *oldgrid, pmegrids_t *newgrid);

#endif
