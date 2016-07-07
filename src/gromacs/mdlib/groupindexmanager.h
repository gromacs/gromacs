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
 * Declares gmx::GroupIndexManager
 *
 * \author Christian Blau <cblau@gwdg.de>
 */

#ifndef GMX_GroupIndexMANAGER_H
#define GMX_GroupIndexMANAGER_H

#include <memory>

#include "groupindex.h"
#include "gromacs/utility/classhelpers.h"


struct gmx_ga2la_t;

namespace gmx
{

/*! \internal \brief
 * Hands out handles to atom group indices and triggers index recalculation for all groups upon domain decomposition if run in parallel.
 *
 * A group index is a collection of local, global and collective indices of atoms.
 *
 */
class GroupIndexManager
{
    public:
        typedef std::shared_ptr<GroupIndex> GroupIndexHandle;
        /*! \brief Construct an group index manager.
         *
         * \param[in] bParallel group index manager needs to know if it manages sets in parallel
         */
        explicit GroupIndexManager(const bool bParallel);
        ~GroupIndexManager();

        /*! \brief Add a new group index to be managed and give back a handle.
         *
         * \param[in] nat The number of atoms in the group index
         * \param[in] ind The indices of the atoms in the group index
         */
        GroupIndexHandle add(const int number_of_atoms, const int *index);

        /*! \brief Erase and destroy all group indexs that only the manager references.
         *
         * Atoms will still be part of the simulation, only their indices will no longer be managed as group index.
         */
        void clean();

        /*! \brief Get the number of managed group indices. */
        size_t numberOfManagedGroupIndices();

        /*! \brief Trigger recalculation of local and collective indices from ga2la if run in parallel.
         *
         * \throws exception if not run in parallel.
         */
        void setIndicesInDomainDecomposition(const gmx_ga2la_t  *ga2la);

    private:
        class Impl;
        PrivateImplPointer<Impl> impl_;
};



} // namespace gmx

#endif
