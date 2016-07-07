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
 * \internal \brief
 * Declares gmx::GroupIndex
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#ifndef GMX_GroupIndex_H
#define GMX_GroupIndex_H

#include <vector>
#include <memory>
#include "gromacs/utility/classhelpers.h"

struct gmx_ga2la_t;
class GroupIndexManager;

namespace gmx
{

/*! \internal \brief
 * A set of atoms is a collection of local, global and collective indices of atoms with a name.
 */
class GroupIndex
{
    public:

        /*! \relates GroupIndexManager
         * Only GroupIndexManger may create, initialize and destroy group indexs and
         * trigger index updates.
         */
        friend class GroupIndexManager;
        ~GroupIndex();

        /*! \brief The collective index maps indices on this node to the global index.
         *
         * \returns the collective index.
         *
         */
        const std::vector<int> &collectiveIndex() const;

        /*! \brief Yields The global as ,e.g., given in an index file.
         *
         * \returns the global index.
         *
         */
        const std::vector<int> &globalIndex() const;

        /*! \brief The local atom indices.
         *
         * \returns the local index.
         *
         */
        const std::vector<int> &localIndex() const;

        /*! \brief The number of atoms from this group index on this node.
         */
        int numAtomsLocal() const;

        /*! \brief The number of all atoms from this group index on all nodes together.
         */
        int numAtomsGlobal() const;
    private:
        /*! \brief Private constructor so Only GroupIndexManager may construct group indexs.
         */
        GroupIndex();


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
