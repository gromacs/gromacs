/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \file
 * \libinternal \brief
 * Declares gmx::LocalAtomSet.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_domdec
 */
#ifndef GMX_LOCALATOMSET_H
#define GMX_LOCALATOMSET_H

#include <vector>

namespace gmx
{

namespace internal
{
class LocalAtomSetData;
}   // namespace internal
/*! \libinternal \brief
 * A local atom set is a collection of local, global and collective indices of atoms.
 *
 *  Local atom sets are constructed by LocalAtomSetManager.
 */
class LocalAtomSet
{
    public:

        friend class LocalAtomSetManager;

        /*! \brief Maps indices on node (0..num_atoms_local_) to global atom indicices.
         *
         * \returns the collective index.
         */
        const std::vector<int> &collectiveIndex() const;

        /*! \brief
         * Global indices of the atoms in this set.
         *
         * \returns the global index.
         */
        const std::vector<int> &globalIndex() const;

        /*! \brief
         * Local indices of the atoms.
         *
         * Access,e.g., the i-th local atom coordinate of this set by x[local_index_[i]].
         * Constructed and updated every domain-decomposition step
         *
         * \returns the local index.
         */
        const std::vector<int> &localIndex() const;
        /*! \brief The number of atoms from this group index on this node.*/
        std::size_t numAtomsLocal() const;
        /*! \brief The number of all atoms from this group index on all nodes together. */
        std::size_t numAtomsGlobal() const;
    private:
        /*! \brief Contructs a new atom set by setting a reference to its internal data.
         * \param[in] data The data for the atom set is stored in a seperate object.
         */
        explicit LocalAtomSet(const internal::LocalAtomSetData &data);

        const internal::LocalAtomSetData * data_;

};

} // namespace gmx

#endif
