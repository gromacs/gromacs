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
// Created by Eric Irrgang on 9/23/17.
//

#include "restraintpotential.h"

#include <cassert>

#include <vector>

#include "gromacs/compat/make_unique.h"

class PotentialContainer::Impl
{
    public:
        std::vector < std::shared_ptr < gmx::IRestraintPotential>> pullers_;
};

void PotentialContainer::addPotential(std::shared_ptr<gmx::IRestraintPotential> puller) noexcept
{
    assert(impl_ != nullptr);
    impl_->pullers_.emplace_back(std::move(puller));
}

PotentialContainer::RestraintIterator PotentialContainer::begin()
{
    return impl_->pullers_.begin();
}

PotentialContainer::RestraintIterator PotentialContainer::end()
{
    return impl_->pullers_.end();
}

PotentialContainer::PotentialContainer() :
    impl_ {gmx::compat::make_unique<PotentialContainer::Impl>()}
{}

//template<typename T>
//std::function<gmx::PotentialPointData(const gmx::Vector &,
//                                           const gmx::Vector &,
//                                           gmx::Time)> gmx::RestraintPotential<T>::getEvaluator()
//{
//    return nullptr;
//}

PotentialContainer::~PotentialContainer() = default;

PotentialContainer &PotentialContainer::operator=(PotentialContainer &&) noexcept = default;

PotentialContainer::PotentialContainer(PotentialContainer &&) noexcept = default;

namespace gmx
{

LegacyPuller::LegacyPuller(pull_t *pullWorkPointer) : pullWorkPointer_(pullWorkPointer)
{};

LegacyPuller::LegacyPuller(const LegacyPuller &source) : LegacyPuller(source.pullWorkPointer_)
{

}

LegacyPuller &LegacyPuller::operator=(const LegacyPuller &source)
{
    this->pullWorkPointer_ = source.pullWorkPointer_;
    return *this;
}

LegacyPuller::LegacyPuller(LegacyPuller &&old) noexcept : LegacyPuller(old.pullWorkPointer_)
{
    old.pullWorkPointer_ = nullptr;
}

LegacyPuller &LegacyPuller::operator=(LegacyPuller &&old) noexcept
{
    this->pullWorkPointer_ = old.pullWorkPointer_;
    old.pullWorkPointer_   = nullptr;
    return *this;
}

PotentialPointData LegacyPuller::calculate(Vector r1,
                                           Vector r2,
                                           Time   t)
{
    (void)(r1);
    (void)(r2);
    (void)(t);
    return {};
}

struct pull_t *LegacyPuller::getRaw()
{
    return pullWorkPointer_;
}

const struct pull_t *LegacyPuller::getRaw() const
{
    return pullWorkPointer_;
}
} // end namespace gmx
