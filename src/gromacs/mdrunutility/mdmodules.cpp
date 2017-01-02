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
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

class MDModules::Impl : public IMDOutputProvider
{
    public:

        Impl()
            : field_(createElectricFieldModule())
        {
        }

        void makeModuleOptions(Options *options)
        {
            // Create a section for applied-forces modules
            auto appliedForcesOptions = options->addSection(OptionSection("applied-forces"));
            field_->mdpOptionProvider()->initMdpOptions(&appliedForcesOptions);
            // In future, other sections would also go here.
        }

        // From IMDOutputProvider
        virtual void initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                                bool bAppendFiles, const gmx_output_env_t *oenv)
        {
            field_->outputProvider()->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);
        }
        virtual void finishOutput()
        {
            field_->outputProvider()->finishOutput();
        }

        std::unique_ptr<IMDModule>      field_;
        std::unique_ptr<ForceProviders> forceProviders_;
};

MDModules::MDModules() : impl_(new Impl)
{
}

MDModules::~MDModules()
{
}

void MDModules::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    auto appliedForcesScope = rules->scopedTransform("/applied-forces");
    impl_->field_->mdpOptionProvider()->initMdpTransform(appliedForcesScope.rules());
}

void MDModules::buildMdpOutput(KeyValueTreeObjectBuilder *builder)
{
    impl_->field_->mdpOptionProvider()->buildMdpOutput(builder);
}

void MDModules::assignOptionsToModules(const KeyValueTreeObject  &params,
                                       IKeyValueTreeErrorHandler *errorHandler)
{
    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);
    // The actual output is in the data fields of the modules that
    // were set up in the module options.
    assignOptionsFromKeyValueTree(&moduleOptions, params, errorHandler);
}

void MDModules::adjustInputrecBasedOnModules(t_inputrec *ir)
{
    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);

    checkForUnknownOptionsInKeyValueTree(*ir->params, moduleOptions);

    std::unique_ptr<KeyValueTreeObject> params(
            new KeyValueTreeObject(
                    adjustKeyValueTreeFromOptions(*ir->params, moduleOptions)));
    delete ir->params;
    ir->params = params.release();
}

IMDOutputProvider *MDModules::outputProvider()
{
    return impl_.get();
}

ForceProviders *MDModules::initForceProviders()
{
    GMX_RELEASE_ASSERT(impl_->forceProviders_ == nullptr,
                       "Force providers initialized multiple times");
    impl_->forceProviders_.reset(new ForceProviders);
    impl_->field_->initForceProviders(impl_->forceProviders_.get());
    return impl_->forceProviders_.get();
}

} // namespace gmx
