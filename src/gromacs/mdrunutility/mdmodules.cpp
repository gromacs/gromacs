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
#include "gmxpre.h"

#include "mdmodules.h"

#include <memory>
#include <vector>

#include "gromacs/applied-forces/electricfield.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

//! Convenience typedef.
using IInputRecExtensionPtr = std::unique_ptr<IInputRecExtension>;

class MDModules::Impl
{
    public:

        Impl() : ir_(nullptr)
        {
        }
        ~Impl()
        {
            if (ir_ != nullptr)
            {
                done_inputrec(ir_);
                sfree(ir_);
            }
        }

        // TODO The set of methods that have to call this do so, but
        // it is not very clear which should do so, rather than assert
        // that someone else should have called this. Perhaps this
        // could be made more intuitive, somehow. Do we need a factory
        // function for MDModules?
        void ensureModulesCreated()
        {
            if (!modulesCreated_)
            {
                modules_.emplace_back(createElectricFieldModule());
                // TODO Give ElectricField an Impl class so the header
                // can reflect the inheritance and this cast is not
                // required.
                // TODO Should we (instead?) have a mechanism where
                // IInputRecExtension objects can register themselves
                // when they are IForceProviders?
                forceProviders_.push_back(dynamic_cast<IForceProvider *>(modules_.back().get()));
                modulesCreated_ = true;
            }
        }

        void ensureInputrecInitialized()
        {
            if (ir_ == nullptr)
            {
                snew(ir_, 1);
                snew(ir_->fepvals, 1);
                snew(ir_->expandedvals, 1);
                snew(ir_->simtempvals, 1);
            }
        }

        bool                               modulesCreated_;
        // TODO Eventually rename IInputRecExtension to IMDModule?
        std::vector<IInputRecExtensionPtr> modules_;
        std::vector<IForceProvider *>      forceProviders_;
        t_inputrec                        *ir_;
        gmx::KeyValueTreeObject            moduleOptionKeyValueTree_;
};

MDModules::MDModules() : impl_(new Impl)
{
}

MDModules::~MDModules()
{
}

t_inputrec *MDModules::inputrec()
{
    impl_->ensureInputrecInitialized();
    return impl_->ir_;
}

const t_inputrec *MDModules::inputrec() const
{
    GMX_RELEASE_ASSERT(impl_->ir_, "Can't access const t_inputrec before creation");
    return impl_->ir_;
}

void MDModules::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    impl_->ensureModulesCreated();
    for(auto &module : impl_->modules_)
    {
        module->initMdpTransform(rules);
    }
}

void MDModules::assignOptionsToModules(KeyValueTreeObject &&optionValues,
                                       IKeyValueTreeErrorHandler *errorHandler)
{
    GMX_RELEASE_ASSERT(impl_->modulesCreated_, "Can't assign options to modules before they have been created");
    impl_->moduleOptionKeyValueTree_ = optionValues;
    inputrec()->moduleOptionParameters = &impl_->moduleOptionKeyValueTree_;
    gmx::Options options;
    for(auto &module : impl_->modules_)
    {
        module->initMdpOptions(&options);
    }
    // TODO Error handling
    gmx::assignOptionsFromKeyValueTree(&options, impl_->moduleOptionKeyValueTree_, errorHandler);
}

void MDModules::printParameters(FILE *fp, int indent) const
{
    GMX_RELEASE_ASSERT(impl_->modulesCreated_, "Can't print parameters unless modules have been created");
    for(auto &module : impl_->modules_)
    {
        module->printParameters(fp, indent);
    }
}

void MDModules::initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                           bool bAppendFiles, const gmx_output_env_t *oenv)
{
    GMX_RELEASE_ASSERT(impl_->modulesCreated_, "Can't start output unless modules have been created");
    for(auto &module : impl_->modules_)
    {
        module->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);
    }
}

void MDModules::finishOutput()
{
    GMX_RELEASE_ASSERT(impl_->modulesCreated_, "Can't finish output unless modules have been created");
    for(auto &module : impl_->modules_)
    {
        module->finishOutput();
    }
}

void MDModules::compare(FILE *fp,
                        const MDModules *other,
                        real reltol,
                        real abstol) const
{
    GMX_RELEASE_ASSERT(impl_->modulesCreated_, "Can't compare modules before they have been created");
    for (auto myIt = std::begin(impl_->modules_), otherIt = std::begin(other->impl_->modules_);
         myIt != std::end(impl_->modules_); ++myIt)
    {
        (*myIt)->compare(fp, otherIt->get(), reltol, abstol);
    }
}

void MDModules::broadCast(const t_commrec *cr)
{
    impl_->ensureModulesCreated();
    for(auto &module : impl_->modules_)
    {
        module->broadCast(cr);
    }
}

void MDModules::initForcerec(t_forcerec *fr)
{
    impl_->ensureModulesCreated();
    for(auto &provider : impl_->forceProviders_)
    {
        provider->initForcerec(fr);
    }
}

void MDModules::calculateForces(const t_commrec  *cr,
                                const t_mdatoms  *mdatoms,
                                PaddedRVecVector *force,
                                double            t)
{
    GMX_RELEASE_ASSERT(impl_->modulesCreated_, "Can't calculate forces for modules before they have been created");
    for(auto &provider : impl_->forceProviders_)
    {
        provider->calculateForces(cr, mdatoms, force, t);
    }
}

} // namespace gmx
