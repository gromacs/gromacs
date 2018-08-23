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
#include <cassert>

#include <iostream>
#include <memory>

#include "md-impl.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"
#include "gmxapi/md/mdmodule.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/keyvaluetree.h"

namespace gmxapi
{

class MDWorkSpec::Impl
{
    public:
        static std::unique_ptr<Impl> create();

        std::vector < std::shared_ptr < gmxapi::MDModule>> modules {};
};

std::unique_ptr<MDWorkSpec::Impl> MDWorkSpec::Impl::create()
{
    auto newImpl = gmx::compat::make_unique<MDWorkSpec::Impl>();
    assert(newImpl != nullptr);
    assert(newImpl->modules.empty());
    return newImpl;
}

MDWorkSpec::MDWorkSpec() :
    impl_ {Impl::create()}
{
    assert(impl_ != nullptr);
}

void MDWorkSpec::addModule(std::shared_ptr<gmxapi::MDModule> module)
{
    assert(impl_ != nullptr);
    std::cout << "Adding module " << module->name() << " to work specification" << std::endl;
    impl_->modules.emplace_back(std::move(module));
}

std::vector < std::shared_ptr < gmxapi::MDModule>> &MDWorkSpec::getModules()
{
    assert(impl_ != nullptr);
    return impl_->modules;
}

MDWorkSpec::~MDWorkSpec() = default;

std::shared_ptr<::gmxapi::MDWorkSpec> MDHolder::getSpec()
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
    return impl_->spec_;
}

std::shared_ptr<const ::gmxapi::MDWorkSpec> MDHolder::getSpec() const
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
    return impl_->spec_;
}

MDHolder::MDHolder() :
    MDHolder {std::make_shared<MDWorkSpec>()}
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
}

MDHolder::MDHolder(std::shared_ptr<MDWorkSpec> spec) :
    name_ {},
impl_ {
    std::make_shared<MDHolder::Impl>(std::move(spec))
}
{
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
}

MDHolder::Impl::Impl(std::shared_ptr<MDWorkSpec> &&spec) :
    spec_ {spec}
{
    assert(spec_ != nullptr);
}

MDHolder::MDHolder(std::string name) :
    MDHolder {}
{
    name_ = name;
    assert(impl_ != nullptr);
    assert(impl_->spec_ != nullptr);
}

std::string MDHolder::name() const
{
    return name_;
}


} //end namespace gmxapi
