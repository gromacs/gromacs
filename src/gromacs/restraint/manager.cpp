/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Implement the restraint manager.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_restraint
 */

#include "gmxpre.h"

#include "manager.h"

#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \internal
 * \brief Implementation class for restraint manager.
 */
class RestraintManager::Impl
{
public:
    /*!
     * \brief Implement Manager::addToSpec()
     *
     * \param restraint Handle to be added to the manager.
     * \param name Identifying string for restraint.
     */
    void add(std::shared_ptr<::gmx::IRestraintPotential> restraint, const std::string& name);

    /*!
     * \brief Clear registered restraints and reset the manager.
     */
    void clear() noexcept;

    /*!
     * \brief The list of configured restraints.
     *
     * Clients can extend the life of a restraint implementation object that
     * is being used by holding the shared_ptr handle. A RestraintManager::Impl is logically
     * const when its owning Manager is logically const, but the Manager can still
     * grant access to individual restraints.
     */
    std::vector<std::shared_ptr<::gmx::IRestraintPotential>> restraint_;

private:
    //! Regulate initialization of the shared resource when (re)initialized.
    static std::mutex initializationMutex_;
};

// Initialize static members
std::mutex RestraintManager::Impl::initializationMutex_{};


void RestraintManager::Impl::add(std::shared_ptr<::gmx::IRestraintPotential> restraint,
                                 const std::string&                          name)
{
    (void)name;
    restraint_.emplace_back(std::move(restraint));
}

RestraintManager::RestraintManager() : instance_(std::make_shared<RestraintManager::Impl>()){};

RestraintManager::~RestraintManager() = default;

void RestraintManager::Impl::clear() noexcept
{
    std::lock_guard<std::mutex> lock(initializationMutex_);
    restraint_.resize(0);
}

void RestraintManager::clear() noexcept
{
    GMX_ASSERT(instance_, "instance_ member should never be null.");
    instance_->clear();
}

void RestraintManager::addToSpec(std::shared_ptr<gmx::IRestraintPotential> puller, const std::string& name)
{
    instance_->add(std::move(puller), name);
}

std::vector<std::shared_ptr<IRestraintPotential>> RestraintManager::getRestraints() const
{
    return instance_->restraint_;
}

unsigned long RestraintManager::countRestraints() noexcept
{
    return instance_->restraint_.size();
}

} // end namespace gmx
