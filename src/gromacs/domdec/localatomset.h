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

#include <memory>
#include <vector>

#include "gromacs/utility/classhelpers.h"

struct gmx_ga2la_t;
class LocalAtomSetManager;

namespace gmx
{

/*! \libinternal \brief
 * A local atom set is a collection of local, global and collective indices of atoms.
 *
 *  Local atom sets are constructed, accessed and destroyed by LocalAtomSetManager.
 */
class LocalAtomSet
{
    public:

        /*! \relates LocalAtomSetManager
         * Only LocalAtomSetManger may create, initialize and destroy group indexs and
         * trigger index updates.
         */
        friend class LocalAtomSetManager;
        ~LocalAtomSet();

        /*! \brief
         * Maps indices on node (0..num_atoms_local_) to global atom indicices.
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
        /*! The number of atoms from this group index on this node.*/
        int numAtomsLocal() const;
        /*! The number of all atoms from this group index on all nodes together. */
        int numAtomsGlobal() const;
    private:

        /*! Private constructor so Only LocalAtomSetManager may construct group indexs. */
        LocalAtomSet();

        /*! \brief Initialize an group index with an index group.
         *
         * \param[in] number_of_atoms The number of atoms of this group index, corresponding to the size of the index
         * \param[in] index An integer index of the atoms in the index set.
         * \param[in] bParallel
         */
        void init(const int number_of_atoms, const int *index, bool bParallel);

        /*! \brief Sets the local and collective indices from a lookup in ga2la.
         *
         * This makes only sense when atoms are indeed distributed over nodes,
         * otherwise local indices equal global indices and are set in init.
         *
         * \param[in] ga2la lookup table that reports if an atom is local.
         */
        void bparallelSetLocalAndCollectiveIndices(const gmx_ga2la_t  *ga2la);

        class Impl;
        PrivateImplPointer<Impl> impl_;

};

} // namespace gmx

#endif
