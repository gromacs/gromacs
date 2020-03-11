/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
    /*! \brief Constructor */
    WholeMoleculeTransform(const gmx_mtop_t& mtop, PbcType pbcType);

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
    //! Buffer for storing coordinates for whole molecules
    std::vector<RVec> wholeMoleculeCoordinates_;
};

} // namespace gmx

#endif
