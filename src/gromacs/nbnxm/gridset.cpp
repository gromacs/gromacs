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
#include "gromacs/utility/fatalerror.h"

#include "atomdata.h"

namespace Nbnxm
{

//! Returns the number of search grids
static int numGrids(const GridSet::DomainSetup& domainSetup)
{
    int numGrids;
    if (domainSetup.doTestParticleInsertion)
    {
        numGrids = 2;
    }
    else
    {
        numGrids = 1;
        for (auto haveDD : domainSetup.haveMultipleDomainsPerDim)
        {
            if (haveDD)
            {
                numGrids *= 2;
            }
        }
    }

    return numGrids;
}

GridSet::DomainSetup::DomainSetup(const int                 ePBC,
                                  const bool                doTestParticleInsertion,
                                  const ivec*               numDDCells,
                                  const gmx_domdec_zones_t* ddZones) :
    ePBC(ePBC),
    doTestParticleInsertion(doTestParticleInsertion),
    haveMultipleDomains(numDDCells != nullptr),
    zones(ddZones)
{
    for (int d = 0; d < DIM; d++)
    {
        haveMultipleDomainsPerDim[d] = (numDDCells != nullptr && (*numDDCells)[d] > 1);
    }
}

GridSet::GridSet(const int                 ePBC,
                 const bool                doTestParticleInsertion,
                 const ivec*               numDDCells,
                 const gmx_domdec_zones_t* ddZones,
                 const PairlistType        pairlistType,
                 const bool                haveFep,
                 const int                 numThreads,
                 gmx::PinningPolicy        pinningPolicy) :
    domainSetup_(ePBC, doTestParticleInsertion, numDDCells, ddZones),
    grids_(numGrids(domainSetup_), Grid(pairlistType, haveFep_)),
    haveFep_(haveFep),
    numRealAtomsLocal_(0),
    numRealAtomsTotal_(0),
    gridWork_(numThreads)
{
    clear_mat(box_);
    changePinningPolicy(&gridSetData_.cells, pinningPolicy);
    changePinningPolicy(&gridSetData_.atomIndices, pinningPolicy);
}

void GridSet::setLocalAtomOrder()
{
    /* Set the atom order for the home cell (index 0) */
    const Nbnxm::Grid& grid = grids_[0];

    int atomIndex = 0;
    for (int cxy = 0; cxy < grid.numColumns(); cxy++)
    {
        const int numAtoms  = grid.numAtomsInColumn(cxy);
        int       cellIndex = grid.firstCellInColumn(cxy) * grid.geometry().numAtomsPerCell;
        for (int i = 0; i < numAtoms; i++)
        {
            gridSetData_.atomIndices[cellIndex] = atomIndex;
            gridSetData_.cells[atomIndex]       = cellIndex;
            atomIndex++;
            cellIndex++;
        }
    }
}

void GridSet::putOnGrid(const matrix                   box,
                        const int                      gridIndex,
                        const rvec                     lowerCorner,
                        const rvec                     upperCorner,
                        const gmx::UpdateGroupsCog*    updateGroupsCog,
                        const gmx::Range<int>          atomRange,
                        real                           atomDensity,
                        gmx::ArrayRef<const int>       atomInfo,
                        gmx::ArrayRef<const gmx::RVec> x,
                        const int                      numAtomsMoved,
                        const int*                     move,
                        nbnxn_atomdata_t*              nbat)
{
    Nbnxm::Grid& grid = grids_[gridIndex];

    int cellOffset;
    if (gridIndex == 0)
    {
        cellOffset = 0;
    }
    else
    {
        const Nbnxm::Grid& previousGrid = grids_[gridIndex - 1];
        cellOffset = previousGrid.atomIndexEnd() / previousGrid.geometry().numAtomsPerCell;
    }

    const int n = atomRange.size();

    real maxAtomGroupRadius;
    if (gridIndex == 0)
    {
        copy_mat(box, box_);

        numRealAtomsLocal_ = *atomRange.end() - numAtomsMoved;
        /* We assume that nbnxn_put_on_grid is called first
         * for the local atoms (gridIndex=0).
         */
        numRealAtomsTotal_ = *atomRange.end() - numAtomsMoved;

        maxAtomGroupRadius = (updateGroupsCog ? updateGroupsCog->maxUpdateGroupRadius() : 0);

        if (debug)
        {
            fprintf(debug, "natoms_local = %5d atom_density = %5.1f\n", numRealAtomsLocal_, atomDensity);
        }
    }
    else
    {
        const Nbnxm::Grid::Dimensions& dimsGrid0 = grids_[0].dimensions();
        atomDensity                              = dimsGrid0.atomDensity;
        maxAtomGroupRadius                       = dimsGrid0.maxAtomGroupRadius;

        numRealAtomsTotal_ = std::max(numRealAtomsTotal_, *atomRange.end());
    }

    /* We always use the home zone (grid[0]) for setting the cell size,
     * since determining densities for non-local zones is difficult.
     */
    const int ddZone = (domainSetup_.doTestParticleInsertion ? 0 : gridIndex);
    // grid data used in GPU transfers inherits the gridset pinning policy
    auto pinPolicy = gridSetData_.cells.get_allocator().pinningPolicy();
    grid.setDimensions(ddZone, n - numAtomsMoved, lowerCorner, upperCorner, atomDensity,
                       maxAtomGroupRadius, haveFep_, pinPolicy);

    for (GridWork& work : gridWork_)
    {
        work.numAtomsPerColumn.resize(grid.numColumns() + 1);
    }

    /* Make space for the new cell indices */
    gridSetData_.cells.resize(*atomRange.end());

    const int nthread = gmx_omp_nthreads_get(emntPairsearch);
    GMX_ASSERT(nthread > 0, "We expect the OpenMP thread count to be set");

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            Grid::calcColumnIndices(grid.dimensions(), updateGroupsCog, atomRange, x, ddZone, move, thread,
                                    nthread, gridSetData_.cells, gridWork_[thread].numAtomsPerColumn);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    /* Copy the already computed cell indices to the grid and sort, when needed */
    grid.setCellIndices(ddZone, cellOffset, &gridSetData_, gridWork_, atomRange, atomInfo.data(), x,
                        numAtomsMoved, nbat);

    if (gridIndex == 0)
    {
        nbat->natoms_local = nbat->numAtoms();
    }
    if (gridIndex == gmx::ssize(grids_) - 1)
    {
        /* We are done setting up all grids, we can resize the force buffers */
        nbat->resizeForceBuffers();
    }

    int maxNumColumns = 0;
    for (int i = 0; i <= gridIndex; i++)
    {
        maxNumColumns = std::max(maxNumColumns, grids_[i].numColumns());
    }
    setNumColumnsMax(maxNumColumns);
}

} // namespace Nbnxm
