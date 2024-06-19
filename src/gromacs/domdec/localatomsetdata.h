/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \file
 * \internal \brief
 * Declares gmx::internal::LocalAtomSetData.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_LOCALATOMSETDATA_H
#define GMX_DOMDEC_LOCALATOMSETDATA_H

#include <numeric>
#include <type_traits>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

class gmx_ga2la_t;

namespace gmx
{

namespace internal
{

/* \brief Internal class for storing and managing atom indices of an atom set.
 */
class LocalAtomSetData
{
public:
    /*! \brief Store the data for an atom set with an index group.
     *
     * Prior to domain decomposition, local atom indices are global atom indices
     * and the collective index runs from 0..numberOfAtoms-1.
     * local and collective indices will be updated in setLocalAndCollectiveIndices
     * to match domain decomposition if domain decomposition is performed.
     *
     * \todo remove this constructor once all indices are represented
     *       as gmx::Index instead of int.
     *
     * \note Not created if the internal int type does match gmx::Index
     *
     * \param[in] globalAtomIndex Indices of the atoms to be managed
     */
    template<typename T = void, typename U = std::enable_if_t<!std::is_same_v<int, Index>, T>>
    explicit LocalAtomSetData(ArrayRef<const int> globalAtomIndex) :
        globalIndex_(globalAtomIndex.begin(), globalAtomIndex.end()),
        localIndex_(globalAtomIndex.begin(), globalAtomIndex.end())
    {
        collectiveIndex_.resize(localIndex_.size());
        std::iota(collectiveIndex_.begin(), collectiveIndex_.end(), 0);
    }

    /*! \brief Store the data for an atom set with an index group.
     *
     * Prior to domain decomposition, local atom indices are global atom indices
     * and the collective index runs from 0..numberOfAtoms-1.
     * local and collective indices will be updated in setLocalAndCollectiveIndices
     * to match domain decomposition if domain decomposition is performed.
     *
     * \param[in] globalAtomIndex Indices of the atoms to be managed
     */
    explicit LocalAtomSetData(ArrayRef<const Index> globalAtomIndex);

    /*! \brief Sets the local and collective indices from a lookup in ga2la.
     *
     * Calculate local and collective indices of home atoms, assuming a valid
     * global atom to local atom look-up table.
     *
     * \param[in] ga2la lookup table that reports if an atom is local.
     */
    void setLocalAndCollectiveIndices(const gmx_ga2la_t& ga2la);
    /*! \brief Global indices of the atoms in this set. */
    const std::vector<int> globalIndex_;
    /*! \brief Maps indices on this rank [0..num_atoms_local_) to global atom indicices,
     * so that localIndex[i] identifies the same atom as globalIndex[collectiveIndex[i]].
     *
     * This translation of locally dense atom data to global representation,
     * allows to adresses per-atom properties, e.g., scattering factors,
     * that are stored in a global continuous array for each atom of the atom set.
     */
    std::vector<int> collectiveIndex_;
    /*! \brief Local indices of the atoms.
     * Access the i-th local atom coordinate of this set by x[local_index_[i]].
     * Constructed and updated every domain-decomposition step.
     */
    std::vector<int> localIndex_;
};

} // namespace internal

} // namespace gmx


#endif
