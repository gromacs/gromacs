/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
#include "gmxpre.h"

#include "mdmodules.h"

#include <cstdio>

#include <memory>
#include <utility>
#include <vector>

#include "gromacs/applied_forces/colvars/colvarsMDModule.h"
#include "gromacs/applied_forces/densityfitting/densityfitting.h"
#include "gromacs/applied_forces/electricfield.h"
#include "gromacs/applied_forces/qmmm/qmmm.h"
#include "gromacs/imd/imd.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

namespace gmx
{

class MDModules::Impl : public IMDOutputProvider
{
public:
    Impl() :
        densityFitting_(DensityFittingModuleInfo::create()),
        field_(createElectricFieldModule()),
        imd_(createInteractiveMolecularDynamicsModule()),
        qmmm_(QMMMModuleInfo::create()),
        swapCoordinates_(createSwapCoordinatesModule()),
        colvars_(ColvarsModuleInfo::create())
    {
    }

    void makeModuleOptions(Options* options) const
    {
        // Create a section for applied-forces modules
        auto appliedForcesOptions = options->addSection(OptionSection("applied-forces"));
        field_->mdpOptionProvider()->initMdpOptions(&appliedForcesOptions);
        densityFitting_->mdpOptionProvider()->initMdpOptions(&appliedForcesOptions);
        qmmm_->mdpOptionProvider()->initMdpOptions(&appliedForcesOptions);
        colvars_->mdpOptionProvider()->initMdpOptions(&appliedForcesOptions);
        // In future, other sections would also go here.
    }

    // From IMDOutputProvider
    void initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv) override
    {
        field_->outputProvider()->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);
        densityFitting_->outputProvider()->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);
    }
    void finishOutput() override
    {
        field_->outputProvider()->finishOutput();
        densityFitting_->outputProvider()->finishOutput();
    }

    /*! \brief Manages callbacks and notifies the MD modules.
     *
     * \note The notifier must be constructed before the modules and shall
     *       not be destructed before the modules are destructed.
     */
    MDModulesNotifiers notifiers_;

    std::unique_ptr<IMDModule>      densityFitting_;
    std::unique_ptr<IMDModule>      field_;
    std::unique_ptr<ForceProviders> forceProviders_;
    std::unique_ptr<IMDModule>      imd_;
    std::unique_ptr<IMDModule>      qmmm_;
    std::unique_ptr<IMDModule>      swapCoordinates_;
    std::unique_ptr<IMDModule>      colvars_;

    /*! \brief List of registered MDModules
     *
     * Note that MDModules::Impl owns this container, but it is only used by
     * the MDModules::initForceProviders() function. To be consistent with
     * IMDModule's vision, as indicated by its docs, we should
     * \todo update IMDModule docs to allow nullptr return values
     * \todo check for nullptr returned by IMDModule methods.
     * \todo include field_ in modules_
     */
    std::vector<std::shared_ptr<IMDModule>> modules_;
};

MDModules::MDModules() : impl_(new Impl) {}

MDModules::~MDModules() {}

void MDModules::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    auto appliedForcesScope = rules->scopedTransform("/applied-forces");
    impl_->field_->mdpOptionProvider()->initMdpTransform(appliedForcesScope.rules());
    impl_->densityFitting_->mdpOptionProvider()->initMdpTransform(appliedForcesScope.rules());
    impl_->qmmm_->mdpOptionProvider()->initMdpTransform(appliedForcesScope.rules());
    impl_->colvars_->mdpOptionProvider()->initMdpTransform(appliedForcesScope.rules());
}

void MDModules::buildMdpOutput(KeyValueTreeObjectBuilder* builder)
{
    impl_->field_->mdpOptionProvider()->buildMdpOutput(builder);
    impl_->densityFitting_->mdpOptionProvider()->buildMdpOutput(builder);
    impl_->qmmm_->mdpOptionProvider()->buildMdpOutput(builder);
    impl_->colvars_->mdpOptionProvider()->buildMdpOutput(builder);
}

void MDModules::assignOptionsToModules(const KeyValueTreeObject& params, IKeyValueTreeErrorHandler* errorHandler)
{
    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);
    // The actual output is in the data fields of the modules that
    // were set up in the module options.
    assignOptionsFromKeyValueTree(&moduleOptions, params, errorHandler);
}

void MDModules::adjustInputrecBasedOnModules(t_inputrec* ir)
{
    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);

    checkForUnknownOptionsInKeyValueTree(*ir->params, moduleOptions);

    std::unique_ptr<KeyValueTreeObject> params(
            new KeyValueTreeObject(adjustKeyValueTreeFromOptions(*ir->params, moduleOptions)));
    delete ir->params;
    ir->params = params.release();
}

IMDOutputProvider* MDModules::outputProvider()
{
    return impl_.get();
}

ForceProviders* MDModules::initForceProviders()
{
    GMX_RELEASE_ASSERT(impl_->forceProviders_ == nullptr,
                       "Force providers initialized multiple times");
    impl_->forceProviders_ = std::make_unique<ForceProviders>();
    impl_->field_->initForceProviders(impl_->forceProviders_.get());
    impl_->densityFitting_->initForceProviders(impl_->forceProviders_.get());
    impl_->qmmm_->initForceProviders(impl_->forceProviders_.get());
    impl_->colvars_->initForceProviders(impl_->forceProviders_.get());
    for (auto&& module : impl_->modules_)
    {
        module->initForceProviders(impl_->forceProviders_.get());
    }
    return impl_->forceProviders_.get();
}

void MDModules::subscribeToPreProcessingNotifications()
{
    impl_->densityFitting_->subscribeToPreProcessingNotifications(&impl_->notifiers_);
    impl_->qmmm_->subscribeToPreProcessingNotifications(&impl_->notifiers_);
    impl_->colvars_->subscribeToPreProcessingNotifications(&impl_->notifiers_);
}

void MDModules::subscribeToSimulationSetupNotifications()
{
    impl_->densityFitting_->subscribeToSimulationSetupNotifications(&impl_->notifiers_);
    impl_->qmmm_->subscribeToSimulationSetupNotifications(&impl_->notifiers_);
    impl_->colvars_->subscribeToSimulationSetupNotifications(&impl_->notifiers_);
}

void MDModules::add(std::shared_ptr<gmx::IMDModule> module)
{
    impl_->modules_.emplace_back(std::move(module));
}

const MDModulesNotifiers& MDModules::notifiers()
{
    return impl_->notifiers_;
}

} // namespace gmx
