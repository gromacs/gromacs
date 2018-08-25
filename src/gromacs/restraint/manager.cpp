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
//
// Created by Eric Irrgang on 10/25/17.
//

/*! \internal \file
 * \brief Implement the restraint manager.
 *
 * \ingroup module_restraint
 */

#include "gmxpre.h"

#include "manager.h"

#include <cassert>

#include <iostream>
#include <map>
#include <memory>
#include <mutex>

#include "gromacs/compat/make_unique.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

namespace restraint
{

// Initialize static members
std::shared_ptr<Manager> Manager::instance_ {
    nullptr
};
std::mutex Manager::initializationMutex_ {};


/*! \internal
 * \brief Implementation class for restraint manager.
 */
class ManagerImpl
{
    public:

        /*!
         * \brief Implement Manager::addToSpec()
         *
         * \param restraint Handle to be added to the manager.
         * \param name Identifying string for restraint.
         */
        void add(std::shared_ptr<::gmx::IRestraintPotential> restraint, std::string name);

        /*!
         * \brief The list of configured restraints.
         *
         * Clients can extend the life of a restraint implementation object that
         * is being used by holding the shared_ptr handle. A ManagerImpl is logically
         * const when its owning Manager is logically const, but the Manager can still
         * grant access to individual restraints.
         */
        mutable std::vector < std::shared_ptr < ::gmx::IRestraintPotential>> restraint_;
};

void ManagerImpl::add(std::shared_ptr<::gmx::IRestraintPotential> restraint, std::string name)
{
    (void)name;
    restraint_.emplace_back(std::move(restraint));
}


Manager::Manager() : impl_(gmx::compat::make_unique<ManagerImpl>()) {};

Manager::~Manager() = default;

void Manager::clear() noexcept
{
    auto new_impl = gmx::compat::make_unique<ManagerImpl>();
    if (new_impl)
    {
        impl_.swap(new_impl);
    }
}

std::shared_ptr<Manager> Manager::instance()
{
    std::lock_guard<std::mutex> lock(initializationMutex_);
    if (instance_ == nullptr)
    {
        // What do we want to do if `new` throws?
        instance_ = std::shared_ptr<Manager>(new Manager);
    }
    assert(instance_ != nullptr);
    return instance_;
}


void Manager::addToSpec(std::shared_ptr<gmx::IRestraintPotential> puller,
                        std::string                               name)
{
    assert(impl_ != nullptr);
    impl_->add(std::move(puller), name);
}

std::vector < std::shared_ptr < IRestraintPotential>> Manager::getSpec() const
{
    if (!impl_)
    {
        GMX_THROW(InternalError("You found a bug!!!\n"
                                "gmx::restraint::Manager::getSpec should not be called before initializing manager!"));
    }
    return impl_->restraint_;
}

unsigned long Manager::countRestraints() noexcept
{
    return impl_->restraint_.size();
}

} // end namespace restraint
} // end namespace gmx
