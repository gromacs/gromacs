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
 * Implements classes in atomsetmanager.h.
 *
 * \author Christian Blau <cblau@gwdg.de>
 */

#include "atomsetmanager.h"
#include "atomset.h"

#include "gromacs/utility/exceptions.h"

#include <algorithm>
#include <memory>

namespace gmx
{

/********************************************************************
 * AtomSetManager::Impl
 */

/*! \internal \brief
 * Private implementation class for AtomSetManager.
 */
class AtomSetManager::Impl
{
    public:
        Impl(bool bParallel);
        bool                bParallel_;
        typedef std::string AtomSetName;
        std::map < AtomSetName, std::unique_ptr < AtomSet>> atom_sets_;
};

AtomSetManager::Impl::Impl(bool bParallel) : bParallel_(bParallel)
{};

/********************************************************************
 * AtomSetManager
 */

/*! \internal \brief
 *  implementation class for TrajectoryAnalysisModuleData.
 *
 * \ingroup module_trajectoryanalysis
 */

AtomSetManager::AtomSetManager(bool bParallel) : impl_(new Impl(bParallel))
{};

AtomSetManager::~AtomSetManager(){};

void AtomSetManager::erase(const std::string &atom_set_name)
{
    if (impl_->atom_sets_.count(atom_set_name) == 0)
    {
        GMX_THROW(InternalError("Cannot erase '" + atom_set_name + "'. Not found in atom sets."));
    }
    impl_->atom_sets_.erase(atom_set_name);
};

void AtomSetManager::add(const std::string &atom_set_name, const int number_of_atoms, const int *index)
{
    if (impl_->atom_sets_.count(atom_set_name) != 0)
    {
        GMX_THROW(InternalError("An atom set with this name already exists. Will not silently overwrite atom sets, use erase(atom_set_name) first."));
    }
    impl_->atom_sets_[atom_set_name] = AtomSet::create();
    impl_->atom_sets_[atom_set_name]->init(number_of_atoms, index, impl_->bParallel_);
};

void AtomSetManager::set_indices_in_domain_decomposition(const gmx_ga2la_t  *ga2la)
{
    if (!impl_->bParallel_)
    {
        GMX_THROW(InternalError("Atom set manager may only set indices in domain decomposition when running in parallel."));
    }
    for (const auto &atom_set : impl_->atom_sets_)
    {
        atom_set.second->bparallel_set_local_and_collective_indices(ga2la);
    }
};

const AtomSet &AtomSetManager::get(const std::string &atom_set_name) const
{
    if (impl_->atom_sets_.count(atom_set_name) == 0)
    {
        GMX_THROW(InternalError("Tried to access atom set whose name is unknown to the atom set manager."));
    }
    return *(impl_->atom_sets_.at(atom_set_name));
};

} // namespace gmx
