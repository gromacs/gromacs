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

#include <cstddef>
#include <cstdint>

#include <array>
#include <memory>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"

#include "grid.h"
#include "gridsetdata.h"

enum class PbcType : int;

namespace gmx
{
class DomdecZones;
enum class PairlistType;
class UpdateGroupsCog;
enum class PinningPolicy : int;
struct nbnxn_atomdata_t;

/*! \internal
 * \brief Holds a set of search grids for the local + non-local DD zones
 *
 * The are three different possible setups:
 * - a single grid, this is the standard case without domain decomposition
 * - one grid for each domain decomposition zone
 * - with test particle insertion there are two grids, one for the system
 *   to insert in and one for the molecule that is inserted
 */
class GridSet
{
public:
    /*! \internal
     * \brief Description of the domain setup: PBC and the connections between domains
     */
    struct DomainSetup
    {
        //! Constructor, without DD \p numDDCells and \p ddZones should be nullptr
        DomainSetup(PbcType            pbcType,
                    bool               doTestParticleInsertion,
                    const IVec*        numDDCells,
                    const DomdecZones* ddZones);

        //! The type of PBC
        PbcType pbcType_;
        //! Tells whether we are doing test-particle insertion
        bool doTestParticleInsertion_;
        //! Are there multiple domains?
        bool haveMultipleDomains;
        //! Are there multiple domains along each dimension?
        std::array<bool, DIM> haveMultipleDomainsPerDim;
        //! The domain decomposition zone setup
        const gmx::DomdecZones* zones;
    };

    //! Constructs a grid set for 1 or multiple DD zones, when numDDCells!=nullptr
    GridSet(PbcType            pbcType,
            bool               doTestParticleInsertion,
            const IVec*        numDDCells,
            const DomdecZones* ddZones,
            PairlistType       pairlistType,
            bool               haveFep,
            int                numThreads,
            PinningPolicy      pinningPolicy);

    //! Puts the atoms on the grid with index \p gridIndex and copies the coordinates to \p nbat
    void putOnGrid(const matrix            box,
                   int                     gridIndex,
                   const rvec              lowerCorner,
                   const rvec              upperCorner,
                   const UpdateGroupsCog*  updateGroupsCog,
                   Range<int>              atomRange,
                   int                     numGridAtoms,
                   real                    atomDensity,
                   ArrayRef<const int32_t> atomInfo,
                   ArrayRef<const RVec>    x,
                   const int*              move,
                   nbnxn_atomdata_t*       nbat);

    //! Returns the domain setup
    DomainSetup domainSetup() const { return domainSetup_; }

    //! Returns the total number of atoms in the grid set, including padding
    int numGridAtomsTotal() const { return grids_.back().atomIndexEnd(); }

    //! Returns the number of local real atoms, i.e. without padded atoms
    int numRealAtomsLocal() const { return numRealAtomsLocal_; }

    //! Returns the number of total real atoms, i.e. without padded atoms
    int numRealAtomsTotal() const { return numRealAtomsTotal_; }

    //! Returns the atom order on the grid for the local atoms
    ArrayRef<const int> getLocalAtomorder() const
    {
        /* Return the atom order for the home cell (index 0) */
        const int numIndices = grids_[0].atomIndexEnd() - grids_[0].firstAtomInColumn(0);

        return constArrayRefFromArray(atomIndices().data(), numIndices);
    }

    //! Sets the order of the local atoms to the order grid atom ordering
    void setLocalAtomOrder();

    //! Return a single grid
    const Grid& grid(size_t idx) const { return grids_[idx]; }

    //! Returns the list of grids
    ArrayRef<const Grid> grids() const { return grids_; }

    //! Returns the grid atom indices covering all grids
    ArrayRef<const int> cells() const { return gridSetData_.cells; }

    //! Returns the grid atom indices covering all grids
    ArrayRef<const int> atomIndices() const { return gridSetData_.atomIndices; }

    //! Returns whether we have perturbed non-bonded interactions
    bool haveFep() const { return haveFep_; }

    //! Returns the unit cell in \p box
    void getBox(matrix box) const { copy_mat(box_, box); }

    //! Returns the maximum number of columns across all grids
    int numColumnsMax() const { return numColumnsMax_; }

    //! Sets the maximum number of columns across all grids
    void setNumColumnsMax(int numColumnsMax) { numColumnsMax_ = numColumnsMax; }

private:
    /* Data members */
    //! The domain setup
    DomainSetup domainSetup_;
    //! The search grids
    std::vector<Grid> grids_;
    //! The cell and atom index data which runs over all grids
    GridSetData gridSetData_;
    //! Tells whether we have perturbed non-bonded interactions
    bool haveFep_;
    //! The periodic unit-cell
    matrix box_;
    //! The number of local real atoms, i.e. without padded atoms, local atoms: 0 to numAtomsLocal_
    int numRealAtomsLocal_;
    //! The total number of real atoms, i.e. without padded atoms
    int numRealAtomsTotal_;
    //! Working data for constructing a single grid, one entry per thread
    std::vector<GridWork> gridWork_;
    //! Maximum number of columns across all grids
    int numColumnsMax_;
};

} // namespace gmx

#endif
