/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/* TODO find out what this file should be called */
#ifndef GMX_EWALD_PME_GRID_H
#define GMX_EWALD_PME_GRID_H

#include <tuple>
#include <vector>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_pme_t;
struct PmeAndFftGrids;

namespace gmx
{
template<typename>
class ArrayRef;
}

template<typename T>
using AlignedVector = std::vector<T, gmx::AlignedAllocator<T>>;

/*! \brief
 * We allow coordinates to be out the unit-cell by up to 2 box lengths,
 * which might be needed along dimension x for a very skewed unit-cell.
 */
constexpr int c_pmeMaxUnitcellShift = 2;

/*! \brief
 * This affects the size of the lookup table of the modulo operation result,
 * when working with PME local grid indices of the particles.
 */
constexpr int c_pmeNeighborUnitcellCount = 2 * c_pmeMaxUnitcellShift + 1;

struct pmegrid_t;
struct pmegrids_t;

void gmx_sum_qgrid_dd(gmx_pme_t* pme, gmx::ArrayRef<real> grid, int direction);

int copy_pmegrid_to_fftgrid(const gmx_pme_t* pme, PmeAndFftGrids* grids);

int copy_fftgrid_to_pmegrid(const gmx_pme_t* pme, PmeAndFftGrids* grids, int nthread, int thread);

void wrap_periodic_pmegrid(const gmx_pme_t* pme, gmx::ArrayRef<real> pmegrid);

void unwrap_periodic_pmegrid(gmx_pme_t* pme, gmx::ArrayRef<real> pmegrid);

/*! \brief Initialized a PME grid struct
 *
 * The actual storage for the grids is passed through \p gridsStorage. This should have
 * size 1 when a single thread is used and 1+nthread with OpenMP threading.
 * When the vectors are empty, they are allocated. When they are not empty the are used
 * and sufficient size is asserted upon.
 */
void pmegrids_init(pmegrids_t*                        grids,
                   int                                nx,
                   int                                ny,
                   int                                nz,
                   int                                nz_base,
                   int                                pme_order,
                   gmx_bool                           bUseThreads,
                   int                                nthread,
                   int                                overlap_x,
                   int                                overlap_y,
                   gmx::ArrayRef<AlignedVector<real>> gridsStorage);

std::tuple<std::vector<int>, std::vector<real>>
make_gridindex_to_localindex(int n, int local_start, int local_range, bool checkRoundingAtBoundary);

void set_grid_alignment(int* pmegrid_nz, int pme_order);

#endif
