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
 * \brief
 * Implements the GridSet class.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "gridset.h"

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/utility/fatalerror.h"

namespace Nbnxm
{

//! Returns the number of DD zones
static int numDDZones(const std::array<bool, DIM> &haveMultipleDomainsPerDim)
{
    int  numDDZones = 1;
    for (auto haveDD : haveMultipleDomainsPerDim)
    {
        if (haveDD)
        {
            numDDZones *= 2;
        }
    }

    return numDDZones;
}

GridSet::DomainSetup::DomainSetup(const int                 ePBC,
                                  const ivec               *numDDCells,
                                  const gmx_domdec_zones_t *ddZones) :
    ePBC(ePBC),
    haveMultipleDomains(numDDCells != nullptr),
    zones(ddZones)
{
    for (int d = 0; d < DIM; d++)
    {
        haveMultipleDomainsPerDim[d] = (numDDCells != nullptr && (*numDDCells)[d] > 1);
    }
}

GridSet::GridSet(const int                 ePBC,
                 const ivec               *numDDCells,
                 const gmx_domdec_zones_t *ddZones,
                 const PairlistType        pairlistType,
                 const bool                haveFep,
                 const int                 numThreads) :
    domainSetup_(ePBC, numDDCells, ddZones),
    grids_(numDDZones(domainSetup_.haveMultipleDomainsPerDim), pairlistType),
    haveFep_(haveFep),
    numRealAtomsLocal_(0),
    numRealAtomsTotal_(0),
    gridWork_(numThreads)
{
    clear_mat(box_);
}

void GridSet::setLocalAtomOrder()
{
    /* Set the atom order for the home cell (index 0) */
    const Nbnxm::Grid &grid = grids_[0];

    int                atomIndex = 0;
    for (int cxy = 0; cxy < grid.numColumns(); cxy++)
    {
        const int numAtoms  = grid.numAtomsInColumn(cxy);
        int       cellIndex = grid.firstCellInColumn(cxy)*grid.geometry().numAtomsPerCell;
        for (int i = 0; i < numAtoms; i++)
        {
            atomIndices_[cellIndex] = atomIndex;
            cells_[atomIndex]       = cellIndex;
            atomIndex++;
            cellIndex++;
        }
    }
}

// TODO: Move to gridset.cpp
void GridSet::putOnGrid(const matrix                    box,
                        const int                       ddZone,
                        const rvec                      lowerCorner,
                        const rvec                      upperCorner,
                        const gmx::UpdateGroupsCog     *updateGroupsCog,
                        const int                       atomStart,
                        const int                       atomEnd,
                        real                            atomDensity,
                        gmx::ArrayRef<const int>        atomInfo,
                        gmx::ArrayRef<const gmx::RVec>  x,
                        const int                       numAtomsMoved,
                        const int                      *move,
                        nbnxn_atomdata_t               *nbat)
{
    Nbnxm::Grid  &grid = grids_[ddZone];

    int           cellOffset;
    if (ddZone == 0)
    {
        cellOffset = 0;
    }
    else
    {
        const Nbnxm::Grid &previousGrid = grids_[ddZone - 1];
        cellOffset = previousGrid.atomIndexEnd()/previousGrid.geometry().numAtomsPerCell;
    }

    const int n = atomEnd - atomStart;

    real      maxAtomGroupRadius;
    if (ddZone == 0)
    {
        copy_mat(box, box_);

        numRealAtomsLocal_ = atomEnd - numAtomsMoved;
        /* We assume that nbnxn_put_on_grid is called first
         * for the local atoms (ddZone=0).
         */
        numRealAtomsTotal_ = atomEnd - numAtomsMoved;

        maxAtomGroupRadius = (updateGroupsCog ? updateGroupsCog->maxUpdateGroupRadius() : 0);

        if (debug)
        {
            fprintf(debug, "natoms_local = %5d atom_density = %5.1f\n",
                    numRealAtomsLocal_, atomDensity);
        }
    }
    else
    {
        const Nbnxm::Grid::Dimensions &dimsGrid0 = grids_[0].dimensions();
        atomDensity        = dimsGrid0.atomDensity;
        maxAtomGroupRadius = dimsGrid0.maxAtomGroupRadius;

        numRealAtomsTotal_ = std::max(numRealAtomsTotal_, atomEnd);
    }

    /* We always use the home zone (grid[0]) for setting the cell size,
     * since determining densities for non-local zones is difficult.
     */
    grid.setDimensions(ddZone, n - numAtomsMoved,
                       lowerCorner, upperCorner,
                       atomDensity,
                       maxAtomGroupRadius,
                       haveFep_);

    for (GridWork &work : gridWork_)
    {
        work.numAtomsPerColumn.resize(grid.numColumns() + 1);
    }

    /* Make space for the new cell indices */
    cells_.resize(atomEnd);

    const int nthread = gmx_omp_nthreads_get(emntPairsearch);

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            Grid::calcColumnIndices(grid.dimensions(),
                                    updateGroupsCog,
                                    atomStart, atomEnd, x,
                                    ddZone, move, thread, nthread,
                                    cells_, gridWork_[thread].numAtomsPerColumn);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    GridSetData gridSetData = getGridSetData();

    /* Copy the already computed cell indices to the grid and sort, when needed */
    grid.setCellIndices(ddZone, cellOffset, &gridSetData, gridWork_,
                        atomStart, atomEnd, atomInfo.data(), x, numAtomsMoved, nbat);

    if (ddZone == 0)
    {
        nbat->natoms_local = nbat->numAtoms();
    }
    if (ddZone == gmx::ssize(grids_) - 1)
    {
        /* We are done setting up all grids, we can resize the force buffers */
        nbat->resizeForceBuffers();
    }
}

} // namespace Nbnxm
