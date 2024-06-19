/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief Declares the WholeMolecules class for generating whole molecules
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_WHOLEMOLECULETRANSFORM_H
#define GMX_MDLIB_WHOLEMOLECULETRANSFORM_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/listoflists.h"

class gmx_ga2la_t;
struct gmx_mtop_t;
enum class PbcType : int;

namespace gmx
{

/*! \libinternal
 * \brief This class manages a coordinate buffer with molecules not split
 * over periodic boundary conditions for use in force calculations
 * which require whole molecules.
 *
 * Note: This class should not be used for computation of forces which
 *       have virial contributions through shift forces.
 */
class WholeMoleculeTransform
{
public:
    /*! \brief Constructor
     *
     * \param[in] mtop               The global topology use for getting the connections between atoms
     * \param[in] pbcType            The type of PBC
     * \param[in] useAtomReordering  Whether we will use atom reordering
     */
    WholeMoleculeTransform(const gmx_mtop_t& mtop, PbcType pbcType, bool useAtomReordering);

    /*! \brief Changes the atom order to the one provided
     *
     * This method is called after domain repartitioning.
     * The object should have been constructed with \p useAtomReordering set to \p true.
     *
     * \param[in] globalAtomIndices  The global atom indices for the local atoms, size should be the system size
     * \param[in] ga2la              Global to local atom index lookup (the inverse of \p globalAtomIndices)
     */
    void updateAtomOrder(ArrayRef<const int> globalAtomIndices, const gmx_ga2la_t& ga2la);

    /*! \brief Updates the graph when atoms have been shifted by periodic vectors */
    void updateForAtomPbcJumps(ArrayRef<const RVec> x, const matrix box);

    /*! \brief Create and return coordinates with whole molecules for input coordinates \p x
     *
     * \param[in] x  Input coordinates, should not have periodic displacement compared
     *               with the coordinates passed in the last call to \p updateForAtomPbcJumps().
     * \param[in] box  The current periodic image vectors
     *
     * Note: this operation is not free. If you need whole molecules coordinates
     * more than once during the force calculation, store the result and reuse it.
     */
    ArrayRef<const RVec> wholeMoleculeCoordinates(ArrayRef<const RVec> x, const matrix box);

private:
    //! The type of PBC
    PbcType pbcType_;
    //! The graph
    t_graph graph_;
    //! The atom index at which graphGlobalAtomOrderEdges_ starts
    int globalEdgeAtomBegin_;
    //! The edges for the global atom order
    ListOfLists<int> graphGlobalAtomOrderEdges_;
    //! Buffer for storing coordinates for whole molecules
    std::vector<RVec> wholeMoleculeCoordinates_;
};

} // namespace gmx

#endif
