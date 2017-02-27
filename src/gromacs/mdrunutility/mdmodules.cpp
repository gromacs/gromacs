/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

#include "gromacs/applied-forces/electricfield.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"
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

        Impl() : field_(nullptr), ir_(nullptr)
        {
            snew(ir_, 1);
            snew(ir_->fepvals, 1);
            snew(ir_->expandedvals, 1);
            snew(ir_->simtempvals, 1);
            // TODO Eventually implement a proper IMDModule, to which
            // create*Module() would return a pointer. It might have
            // methods in its interface that return IInputRecExtension
            // (renamed IMdpOptionsProvider) and IForceProvider.
            field_      = createElectricFieldModule();
            ir_->efield = field_.get();
        }
        ~Impl()
        {
            if (ir_ != nullptr)
            {
                done_inputrec(ir_);
                sfree(ir_);
            }
        }

        void makeModuleOptions(Options *options)
        {
            // Create a section for applied-forces modules
            auto appliedForcesOptions = options->addSection(OptionSection("applied-forces"));
            field_->initMdpOptions(&appliedForcesOptions);
            // In future, other sections would also go here.
        }

        IInputRecExtensionPtr  field_;
        t_inputrec            *ir_;
};

MDModules::MDModules() : impl_(new Impl)
{
}

MDModules::~MDModules()
{
}

t_inputrec *MDModules::inputrec()
{
    return impl_->ir_;
}

const t_inputrec *MDModules::inputrec() const
{
    return impl_->ir_;
}

void MDModules::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    // TODO The transform rules for applied-forces modules should
    // embed the necessary prefix (and similarly for other groupings
    // of modules). For now, electric-field embeds this itself.
    impl_->field_->initMdpTransform(rules);
}

void MDModules::assignOptionsToModulesFromMdp(const KeyValueTreeObject  &mdpOptionValues,
                                              IKeyValueTreeErrorHandler *errorHandler)
{
    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);

    KeyValueTreeObject keyValueParameters(mdpOptionValues);
    impl_->ir_->params = new KeyValueTreeObject(adjustKeyValueTreeFromOptions(keyValueParameters, moduleOptions));
    // The actual output is in the data fields of the modules that
    // were set up in the module options.
    assignOptionsFromKeyValueTree(&moduleOptions, *impl_->ir_->params, errorHandler);
}

void MDModules::assignOptionsToModulesFromTpr()
{
    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);

    // Note that impl_->ir_->params was set up during tpr reading, so
    // all we need to do here is integrate that with the module
    // options, which e.g. might have changed between versions.
    // The actual output is in the data fields of the modules that
    // were set up in the module options.
    //
    // TODO error handling
    assignOptionsFromKeyValueTree(&moduleOptions, *impl_->ir_->params, nullptr);
}

void MDModules::adjustInputrecBasedOnModules()
{
    gmx::Options                        options;
    impl_->field_->initMdpOptions(&options);
    std::unique_ptr<KeyValueTreeObject> params(impl_->ir_->params);
    // Avoid double freeing if the next operation throws.
    impl_->ir_->params = nullptr;
    impl_->ir_->params = new KeyValueTreeObject(
                gmx::adjustKeyValueTreeFromOptions(*params, options));
}

} // namespace gmx
