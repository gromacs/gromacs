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
#include <vector>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#include "localatomsetdata.h"

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
    std::vector<std::unique_ptr<internal::LocalAtomSetData>> atomSetData_; /**< handles to the managed atom sets */
};

/********************************************************************
 * LocalAtomSetManager */

LocalAtomSetManager::LocalAtomSetManager() : impl_(new Impl()) {}

LocalAtomSetManager::~LocalAtomSetManager() {}

template<>
LocalAtomSet LocalAtomSetManager::add<void, void>(ArrayRef<const int> globalAtomIndex)
{
    impl_->atomSetData_.push_back(std::make_unique<internal::LocalAtomSetData>(globalAtomIndex));
    return LocalAtomSet(*impl_->atomSetData_.back());
}

LocalAtomSet LocalAtomSetManager::add(ArrayRef<const Index> globalAtomIndex)
{
    impl_->atomSetData_.push_back(std::make_unique<internal::LocalAtomSetData>(globalAtomIndex));
    return LocalAtomSet(*impl_->atomSetData_.back());
}

void LocalAtomSetManager::setIndicesInDomainDecomposition(const gmx_ga2la_t& ga2la)
{
    for (const auto& atomSet : impl_->atomSetData_)
    {
        atomSet->setLocalAndCollectiveIndices(ga2la);
    }
}

} // namespace gmx
