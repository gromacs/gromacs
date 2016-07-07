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
/*! \internal \file
 * \brief
 * Implements classes in LocalAtomSetmanager.h.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "localatomsetmanager.h"

#include <algorithm>
#include <memory>

#include "gromacs/utility/exceptions.h"

#include "localatomset.h"

namespace gmx
{

/********************************************************************
 * LocalAtomSetManager::Impl */

/*! \internal \brief
 * Private implementation class for LocalAtomSetManager.
 */
class LocalAtomSetManager::Impl
{
    public:
        /*! \brief Impl is meaningless without knowing if run in parallel.
         *
         * \param[in] bParallel true if simulation is run in parallel */
        explicit Impl(bool bParallel);
        bool bParallel_;                      /**< set true if atoms are spread over nodes*/
        std::vector<LocalAtomSet> atom_sets_; /**< handles to the managed atom sets */
};

LocalAtomSetManager::Impl::Impl(bool bParallel) : bParallel_(bParallel)
{}

/********************************************************************
 * LocalAtomSetManager */

LocalAtomSetManager::LocalAtomSetManager(bool bParallel) : impl_(new Impl(bParallel))
{}

LocalAtomSetManager::~LocalAtomSetManager(){}

const LocalAtomSet &
LocalAtomSetManager::add(const int number_of_atoms, const int *index)
{
    LocalAtomSet new_set(number_of_atoms, index, impl_->bParallel_);
    impl_->atom_sets_.push_back(std::move(new_set));
    return impl_->atom_sets_.back();
}

void
LocalAtomSetManager::setIndicesInDomainDecomposition(const gmx_ga2la_t  *ga2la)
{
    if (!impl_->bParallel_)
    {
        GMX_THROW(InternalError("group index manager may only set indices in domain decomposition when running in parallel."));
    }
    for (auto &atom_set : impl_->atom_sets_)
    {
        atom_set.setLocalAndCollectiveIndices(ga2la);
    }
}

} // namespace gmx
