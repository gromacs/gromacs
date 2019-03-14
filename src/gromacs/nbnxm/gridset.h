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
 * Declares the GridSet class.
 *
 * This class holds the grids for the local and non-local domain decomposition
 * zones, as well as the cell and atom data that covers all grids.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GRIDSET_H
#define GMX_NBNXM_GRIDSET_H

#include <memory>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"

#include "grid.h"


struct nbnxn_atomdata_t;
enum class PairlistType;

namespace gmx
{
class UpdateGroupsCog;
}

namespace Nbnxm
{

/*! \internal
 * \brief An object holding a set of search grids for the local + non-local DD zones
 */
class GridSet
{
    public:
        //! Constructs a grid set for 1 or multiple DD zones, when numDDCells!=nullptr
        GridSet(const std::array<bool, DIM> &haveDomDecPerDim,
                PairlistType                 pairlistType,
                bool                         haveFep,
                int                          numThreads);

        //! Puts the atoms in \p ddZone on the grid and copies the coordinates to \p nbat
        void putOnGrid(const matrix                    box,
                       int                             ddZone,
                       const rvec                      lowerCorner,
                       const rvec                      upperCorner,
                       const gmx::UpdateGroupsCog     *updateGroupsCog,
                       int                             atomStart,
                       int                             atomEnd,
                       real                            atomDensity,
                       const int                      *atinfo,
                       gmx::ArrayRef<const gmx::RVec>  x,
                       int                             numAtomsMoved,
                       const int                      *move,
                       nbnxn_atomdata_t               *nbat);

        //! Returns the number of cells along x and y for the local grid
        void getLocalNumCells(int *numCellsX,
                              int *numCellsY) const
        {
            *numCellsX = grids_[0].dimensions().numCells[XX];
            *numCellsY = grids_[0].dimensions().numCells[YY];
        }

        //! Returns the total number of atoms in the grid set, including padding
        int numGridAtomsTotal() const
        {
            return grids_.back().atomIndexEnd();
        }

        //! Returns the number of local real atoms, i.e. without padded atoms
        int numRealAtomsLocal() const
        {
            return numRealAtomsLocal_;
        }

        //! Returns the number of total real atoms, i.e. without padded atoms
        int numRealAtomsTotal() const
        {
            return numRealAtomsTotal_;
        }

        //! Returns the atom order on the grid for the local atoms
        gmx::ArrayRef<const int> getLocalAtomorder() const
        {
            /* Return the atom order for the home cell (index 0) */
            const int numIndices = grids_[0].atomIndexEnd() - grids_[0].firstAtomInColumn(0);

            return gmx::constArrayRefFromArray(atomIndices_.data(), numIndices);
        }

        //! Sets the order of the local atoms to the order grid atom ordering
        void setLocalAtomOrder();

        //! Returns the list of grids
        gmx::ArrayRef<const Grid> grids() const
        {
            return grids_;
        }

        //! Returns the grid atom indices covering all grids
        gmx::ArrayRef<const int> cells() const
        {
            return cells_;
        }

        //! Returns the grid atom indices covering all grids
        gmx::ArrayRef<const int> atomIndices() const
        {
            return atomIndices_;
        }

        //! Returns whether we have perturbed non-bonded interactions
        bool haveFep() const
        {
            return haveFep_;
        }

        //! Returns the unit cell in \p box
        void getBox(matrix box) const
        {
            copy_mat(box_, box);
        }

    private:
        //! Returns collection of the data that covers all grids
        const GridSetData getGridSetData()
        {
            GridSetData gridSetData = { cells_, atomIndices_, haveFep_ };

            return gridSetData;
        }

        /* Data members */
        //! The search grids
        std::vector<Grid>     grids_;
        //! The actual cell indices for all atoms, covering all grids
        std::vector<int>      cells_;
        //! The actual array of atom indices, covering all grids
        std::vector<int>      atomIndices_;
        //! Tells whether we have perturbed non-bonded interactions
        bool                  haveFep_;
        //! The periodic unit-cell
        matrix                box_;
        //! The number of local real atoms, i.e. without padded atoms, local atoms: 0 to numAtomsLocal_
        int                   numRealAtomsLocal_;
        //! The total number of real atoms, i.e. without padded atoms
        int                   numRealAtomsTotal_;
        //! Working data for constructing a single grid, one entry per thread
        std::vector<GridWork> gridWork_;
};

} // namespace Nbnxm

#endif
