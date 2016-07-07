/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \inlibraryapi
 * \ingroup module_domdec
 */
#ifndef GMX_LOCALATOMSET_H
#define GMX_LOCALATOMSET_H

#include "localatomsetdata.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{
/*! \libinternal \brief
 * A local atom set collects local, global and collective indices of
 * the home atoms on a rank. The indices of the home atoms are automatically
 * updated during domain decomposition, thus gmx::LocalAtomSet::localIndex
 * enables iteration over local atoms properties like coordinates or forces.
 * TODO: add a LocalAtomSet iterator.
 *
 * To generate a LocalAtomSet call gmx::LocalAtomSetManger::add and keep the
 * handle to the LocalAtomSet returned from this call.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 */
class LocalAtomSet
{
    public:
        friend class LocalAtomSetManager;
        /*! \brief Maps indices on node (0..numAtomsLocal_) to global atom indicices.
         *
         * \returns the collective index.
         */
        ArrayRef<const int> collectiveIndex() const
        {
            return data_.collectiveIndex_;
        }
        /*! \brief Global indices of the atoms in this set.
         *
         * \returns the global index.
         */
        ArrayRef<const int> globalIndex() const
        {
            return data_.globalIndex_;
        }
        /*! \brief Local indices of the atoms.
         *
         * Access,e.g., the i-th local atom coordinate of this set by x[local_index_[i]].
         * Updated within domain-decomposition.
         *
         * \returns the local index.
         */
        ArrayRef<const int> localIndex() const
        {
            return data_.localIndex_;
        }
        /*! \brief The number of atoms from this group index on this node.*/
        std::size_t numAtomsLocal() const
        {
            return data_.localIndex_.size();
        }
        /*! \brief The number of all atoms from this group index on all nodes together. */
        std::size_t numAtomsGlobal() const
        {
            return data_.globalIndex_.size();
        }
    private:
        /*! \brief Contructs a new atom set by setting a reference to its internal data.
         * \param[in] data The data for the atom set is stored
         * in LocalAtomSetData, which is manged by \ref gmx::LocalAtomSetManager.
         */
        explicit LocalAtomSet(const internal::LocalAtomSetData &data) : data_ {data}
        {}
        const internal::LocalAtomSetData &data_;

};

} // namespace gmx

#endif
