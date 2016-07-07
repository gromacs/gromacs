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
 * \brief
 * Declares gmx::AtomSet
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#ifndef GMX_ATOMSET_H
#define GMX_ATOMSET_H

#include <vector>
#include <memory>
#include "gromacs/utility/classhelpers.h"

struct gmx_ga2la_t;
class AtomSetManager;

namespace gmx
{

/*! \internal \brief
 * A set of atoms is a collection of local, global and collective indices of atoms with a name.
 */
class AtomSet
{
    public:


        /*! \brief The collective index maps indices on this node to the global index.
         *
         * \returns the collective index.
         *
         */
        const std::vector<int> &collective_index() const;

        /*! \brief Yields The global as ,e.g., given in an index file.
         *
         * \returns the global index.
         *
         */
        const std::vector<int> &global_index() const;

        /*! \brief The local atom indices.
         *
         * \returns the local index.
         *
         */
        const std::vector<int> &local_index() const;


        /*! \brief The number of atoms from this atom set on this node.
         */
        int num_atoms_local() const;

        /*! \brief The number of all atoms from this atom set on all nodes together.
         */
        int num_atoms_global() const;

        /*! \brief Begin iterator for local index, named begin() so range based for loop over atom set works.
         */
        std::vector<int>::const_iterator begin() const;
        /*! \brief End iterator for local index, named end() so range based for loop over atom set works.
         */
        std::vector<int>::const_iterator end() const;

        /*! \relates AtomSetManager
         * Only AtomSetManger may create, initialize and destroy atom sets and
         * trigger index updates.
         */
        friend class AtomSetManager;

        /* To have a private destructor but still return a unique_ptr reference
         * to AtomSet, befriend deleter_type */
        friend std::unique_ptr<AtomSet>::deleter_type;

    private:
        /*! \brief Private constructor so Only AtomSetManager may construct atom sets.
         * \todo: ideally hide this from AtomSetManager which should only call create().
         */
        AtomSet();

        /*! \brief Private destrcutor so only AtomSetManager may destroy atom sets.
         */
        ~AtomSet();

        /*! \brief Make an unitialized atom set; only AtomSetManager may do that.
         *
         * With this setup we can, e.g., test the atom set manager without having to provide actual indices.
         *
         * \returns A unique pointer, because atom sets are meant to be owned by atom set manager.
         */
        static std::unique_ptr<AtomSet> create();

        /*! \brief Initialize an atom set with an index group.
         *
         * \param[in] number_of_atoms The number of atoms of this atom set, corresponding to the size of the index
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
        void bparallel_set_local_and_collective_indices(const gmx_ga2la_t  *ga2la);

        class Impl;
        PrivateImplPointer<Impl> impl_;

};



} // namespace gmx

#endif
