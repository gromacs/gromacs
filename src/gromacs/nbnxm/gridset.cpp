/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

#include <cmath>
#include <cstdio>

#include <algorithm>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/gridsetdata.h"
#include "gromacs/utility/allocator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

//! Returns the number of search grids
static int numGrids(const GridSet::DomainSetup& domainSetup)
{
    // One grid for the test particle, one for the rest
    static constexpr int sc_numGridsForTestParticleInsertion = 2;
    if (domainSetup.doTestParticleInsertion_)
    {
        return sc_numGridsForTestParticleInsertion;
    }
    else
    {
        int numGrids = 1;
        for (auto haveDD : domainSetup.haveMultipleDomainsPerDim)
        {
            if (haveDD)
            {
                numGrids *= 2;
            }
        }
        return numGrids;
    }
}

GridSet::DomainSetup::DomainSetup(const PbcType      pbcType,
                                  const bool         doTestParticleInsertion,
                                  const IVec*        numDDCells,
                                  const DomdecZones* ddZones) :
    pbcType_(pbcType),
    doTestParticleInsertion_(doTestParticleInsertion),
    haveMultipleDomains(numDDCells != nullptr
                        && (*numDDCells)[XX] * (*numDDCells)[YY] * (*numDDCells)[ZZ] > 1),
    zones(ddZones)
{
    for (int d = 0; d < DIM; d++)
    {
        haveMultipleDomainsPerDim[d] = (numDDCells != nullptr && (*numDDCells)[d] > 1);
    }
}

GridSet::GridSet(const PbcType      pbcType,
                 const bool         doTestParticleInsertion,
                 const IVec*        numDDCells,
                 const DomdecZones* ddZones,
                 const PairlistType pairlistType,
                 const bool         haveFep,
                 const int          numThreads,
                 PinningPolicy      pinningPolicy) :
    domainSetup_(pbcType, doTestParticleInsertion, numDDCells, ddZones),
    grids_(numGrids(domainSetup_), Grid(pairlistType, haveFep_, pinningPolicy)),
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
    const Grid& grid = grids_[0];

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

static int getGridOffset(ArrayRef<const Grid> grids, int gridIndex)
{
    if (gridIndex == 0)
    {
        return 0;
    }
    else
    {
        const Grid& previousGrid = grids[gridIndex - 1];
        return previousGrid.atomIndexEnd() / previousGrid.geometry().numAtomsPerCell;
    }
}

void GridSet::putOnGrid(const matrix            box,
                        const int               gridIndex,
                        const rvec              lowerCorner,
                        const rvec              upperCorner,
                        const UpdateGroupsCog*  updateGroupsCog,
                        const Range<int>        atomRange,
                        const int               numGridAtoms,
                        real                    atomDensity,
                        ArrayRef<const int32_t> atomInfo,
                        ArrayRef<const RVec>    x,
                        const int*              move,
                        nbnxn_atomdata_t*       nbat)
{
    Grid&     grid               = grids_[gridIndex];
    const int cellOffset         = getGridOffset(grids_, gridIndex);
    real      maxAtomGroupRadius = NAN;

    if (gridIndex == 0)
    {
        copy_mat(box, box_);

        numRealAtomsLocal_ = numGridAtoms;
        /* We assume that nbnxn_put_on_grid is called first
         * for the local atoms (gridIndex=0).
         */
        numRealAtomsTotal_ = numGridAtoms;

        maxAtomGroupRadius = (updateGroupsCog ? updateGroupsCog->maxUpdateGroupRadius() : 0);

        if (debug)
        {
            fprintf(debug, "natoms_local = %5d atom_density = %5.1f\n", numRealAtomsLocal_, atomDensity);
        }
    }
    else
    {
        const Grid::Dimensions& dimsGrid0 = grids_[0].dimensions();
        atomDensity                       = dimsGrid0.atomDensity;
        maxAtomGroupRadius                = dimsGrid0.maxAtomGroupRadius;

        numRealAtomsTotal_ = std::max(numRealAtomsTotal_, *atomRange.end());
    }

    /* We always use the home zone (grid[0]) for setting the cell size,
     * since determining densities for non-local zones is difficult.
     */
    const int ddZone = (domainSetup_.doTestParticleInsertion_ ? 0 : gridIndex);

    /* The search grid/cells should be optimized to be as close to cubic
     * as possible. This is not possible to achieve in a single go
     * in cases where the particles are not distributed homogeneously.
     * Thus we put the particles on the 2D-grid in 1 or 2 iterations for zone 0.
     * We first generate a grid that is optimal for a homogeneous particle
     * density. We then compute the effective grid density. If this is more
     * than a factor 1.5 higher than the homogeneous density, we use a finer grid
     * based on the newly computed density.
     */
    const real c_gridDensityRatioThreshold = 1.5_real;
    const bool optimizeDensity             = (ddZone == 0 && !atomRange.empty());
    real       gridDensityRatio            = 0;
    int        iteration                   = 0;

    while (iteration == 0
           || (optimizeDensity && iteration == 1 && gridDensityRatio > c_gridDensityRatioThreshold))
    {
        if (iteration == 1)
        {
            /* The effective 2D grid density is higher than the uniform density.
             * So we need to increase the 3D density, but we only know about
             * the density in 2D. If the inhomogeneity is 2D only (unlikely),
             * we need to correct with a factor gridDensityRatio. If we have
             * a sphere-like concentration of particles, the correction factor
             * should be gridDensityRatio^3/2. We use the average exponent.
             */
            atomDensity *= std::pow(gridDensityRatio, 1.25_real);
        }

        const bool computeGridDensityRatio = (iteration == 0 && optimizeDensity);

        gridDensityRatio = generateAndFill2DGrid(&grid,
                                                 gridWork_,
                                                 &gridSetData_.cells,
                                                 lowerCorner,
                                                 upperCorner,
                                                 updateGroupsCog,
                                                 atomRange,
                                                 numGridAtoms,
                                                 &atomDensity,
                                                 maxAtomGroupRadius,
                                                 x,
                                                 ddZone,
                                                 move,
                                                 computeGridDensityRatio);

        iteration++;
    }

    /* Copy the already computed cell indices to the grid and sort, when needed */
    grid.setCellIndices(
            ddZone, cellOffset, &gridSetData_, gridWork_, atomRange, numGridAtoms, atomInfo, x, nbat);

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

} // namespace gmx
