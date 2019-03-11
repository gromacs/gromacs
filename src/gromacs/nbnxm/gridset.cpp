/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "gromacs/nbnxm/nbnxm.h"

#include "internal.h"

namespace Nbnxm
{

static int numGrids(const ivec *numDDCells)
{
    int  numGrids   = 1;
    bool haveDomDec = (numDDCells != nullptr);
    if (haveDomDec)
    {
        for (int d = 0; d < DIM; d++)
        {
            if ((*numDDCells)[d] > 1)
            {
                /* Each grid matches a DD zone */
                numGrids *= 2;
            }
        }
    }

    return numGrids;
}

GridSet::GridSet(const ivec         *numDDCells,
                 const PairlistType  pairlistType,
                 const bool          haveFep,
                 const int           numThreads) :
    grids_(numGrids(numDDCells), pairlistType),
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

} // namespace Nbnxm


void nonbonded_verlet_t::setLocalAtomOrder()
{
    nbs->gridSet_.setLocalAtomOrder();
}

void nonbonded_verlet_t::getLocalNumCells(int *numCellsX,
                                          int *numCellsY) const
{
    nbs->gridSet().getLocalNumCells(numCellsX, numCellsY);
}

gmx::ArrayRef<const int> nbnxn_get_gridindices(const nbnxn_search* nbs)
{
    return nbs->gridSet().cells();
}
