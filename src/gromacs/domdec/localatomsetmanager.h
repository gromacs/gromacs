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
/*! \libinternal \file
 * \brief
 * Declares gmx::LocalAtomSetManager
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_domdec
 */

#ifndef GMX_LOCALATOMSETMANAGER_H
#define GMX_LOCALATOMSETMANAGER_H

#include <memory>

#include "gromacs/utility/classhelpers.h"

struct gmx_ga2la_t;

namespace gmx
{

class LocalAtomSet;

/*! \libinternal \brief
 * Hands out handles to local atom set indices and triggers index recalculation
 * for all sets upon domain decomposition if run in parallel.
 *
 * \ingroup module_domdec
 */
class LocalAtomSetManager
{
    public:
        /*! \brief Construct an atom set manager.
         *
         * \param[in] bParallel atom set manager needs to know if it manages sets in a parallel run
         */
        explicit LocalAtomSetManager(const bool bParallel);
        ~LocalAtomSetManager();

        /*! \brief Add a new atom set to be managed and give back a handle.
         *
         * \param[in] number_of_atoms The number of atoms in the atom set
         * \param[in] index The indices of the atoms in the atom set
         */
        const LocalAtomSet &add(const int number_of_atoms, const int *index);
        void remove(const LocalAtomSet &atom_set);

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
