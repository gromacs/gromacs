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
#include <unordered_map>
#include <utility>

#include "gromacs/applied_forces/colvars/colvarsMDModule.h"
#include "gromacs/applied_forces/densityfitting/densityfitting.h"
#include "gromacs/applied_forces/electricfield.h"
#include "gromacs/applied_forces/nnpot/nnpot.h"
#include "gromacs/applied_forces/plumed/plumedMDModule.h"
#include "gromacs/applied_forces/qmmm/qmmm.h"
#include "gromacs/fmm/fmm_mdmodule.h"
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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;

namespace gmx
{

class MDModules::Impl : public IMDOutputProvider
{
public:
    Impl()
    {
        modules_[std::string(DensityFittingModuleInfo::sc_name)] = DensityFittingModuleInfo::create();
        modules_[std::string(ElectricFieldModuleInfo::sc_name)] = ElectricFieldModuleInfo::create();
        modules_[std::string(InteractiveMolecularDynamicsModuleInfo::sc_name)] =
                InteractiveMolecularDynamicsModuleInfo::create();
        modules_[std::string(QMMMModuleInfo::sc_name)] = QMMMModuleInfo::create();
        modules_[std::string(SwapCoordinatesModuleInfo::sc_name)] = SwapCoordinatesModuleInfo::create();
        modules_[std::string(ColvarsModuleInfo::sc_name)] = ColvarsModuleInfo::create();
        modules_[std::string(PlumedModuleInfo::sc_name)]  = PlumedModuleInfo::create();
        modules_[std::string(NNPotModuleInfo::sc_name)]   = NNPotModuleInfo::create();
        modules_[std::string(FmmModuleInfo::sc_name)]     = FmmModuleInfo::create();
    }

    void makeModuleOptions(Options* options) const
    {
        // TODO When we have mdp option providers other than
        // applied-forces modules, do we need to extend IMdpOptions to
        // return an optional section name, like "applied-forces"?
        // What value does it provide other than prepending a string?
        // Can we handle it differently?

        // Create a section for applied-forces modules
        auto appliedForcesOptions = options->addSection(OptionSection("applied-forces"));
        for (std::string_view moduleName : { ElectricFieldModuleInfo::sc_name,
                                             DensityFittingModuleInfo::sc_name,
                                             QMMMModuleInfo::sc_name,
                                             ColvarsModuleInfo::sc_name,
                                             NNPotModuleInfo::sc_name })
        {
            IMDModule*          module            = modules_.at(std::string(moduleName)).get();
            IMdpOptionProvider* mdpOptionProvider = module->mdpOptionProvider();
            GMX_RELEASE_ASSERT(mdpOptionProvider,
                               "Applied-forces modules all implement MDP options");
            mdpOptionProvider->initMdpOptions(&appliedForcesOptions);
        }

        auto       fmmModuleOptions = options->addSection(OptionSection("fast-multipole-method"));
        IMDModule* fmmModule        = modules_.at(std::string(FmmModuleInfo::sc_name)).get();
        IMdpOptionProvider* fmmMdpOptionProvider = fmmModule->mdpOptionProvider();
        GMX_RELEASE_ASSERT(fmmMdpOptionProvider,
                           "Fast-multipole-method module must implement MDP options");
        fmmMdpOptionProvider->initMdpOptions(&fmmModuleOptions);

        // In future, other sections would also go here.
    }

    // From IMDOutputProvider
    void initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv) override
    {
        for (auto& [_, module] : modules_)
        {
            if (module->outputProvider())
            {
                module->outputProvider()->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);
            }
        }
    }
    void finishOutput() override
    {
        for (auto& [_, module] : modules_)
        {
            if (module->outputProvider())
            {
                module->outputProvider()->finishOutput();
            }
        }
    }

    /*! \brief Manages callbacks and notifies the MD modules.
     *
     * \note The notifier must be constructed before the modules and shall
     *       not be destructed before the modules are destructed.
     */
    MDModulesNotifiers notifiers_;

    std::unique_ptr<ForceProviders> forceProviders_;

    /*! \brief List of registered MDModules */
    std::unordered_map<std::string, std::shared_ptr<IMDModule>> modules_;
};

MDModules::MDModules() : impl_(new Impl) {}

MDModules::~MDModules() {}

void MDModules::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    auto appliedForcesScope = rules->scopedTransform("/applied-forces");
    for (std::string_view moduleName : { ElectricFieldModuleInfo::sc_name,
                                         DensityFittingModuleInfo::sc_name,
                                         QMMMModuleInfo::sc_name,
                                         ColvarsModuleInfo::sc_name,
                                         NNPotModuleInfo::sc_name })
    {
        IMDModule*          module            = impl_->modules_.at(std::string(moduleName)).get();
        IMdpOptionProvider* mdpOptionProvider = module->mdpOptionProvider();
        GMX_RELEASE_ASSERT(mdpOptionProvider, "Applied-forces modules all implement MDP options");
        mdpOptionProvider->initMdpTransform(appliedForcesScope.rules());
    }

    auto                fmmScope  = rules->scopedTransform("/fast-multipole-method");
    IMDModule*          fmmModule = impl_->modules_.at(std::string(FmmModuleInfo::sc_name)).get();
    IMdpOptionProvider* fmmMdpOptionProvider = fmmModule->mdpOptionProvider();
    GMX_RELEASE_ASSERT(fmmMdpOptionProvider,
                       "Fast-multipole-method module must implement MDP options");
    fmmMdpOptionProvider->initMdpTransform(fmmScope.rules());
}

void MDModules::buildMdpOutput(KeyValueTreeObjectBuilder* builder)
{
    // Note that the order here should be kept stable so that the
    // fields in the mdp output file appear in a reliable order.
    for (std::string_view moduleName : { ElectricFieldModuleInfo::sc_name,
                                         DensityFittingModuleInfo::sc_name,
                                         QMMMModuleInfo::sc_name,
                                         ColvarsModuleInfo::sc_name,
                                         NNPotModuleInfo::sc_name })
    {
        IMDModule*                module = impl_->modules_.at(std::string(moduleName)).get();
        const IMdpOptionProvider* mdpOptionProvider = module->mdpOptionProvider();
        GMX_RELEASE_ASSERT(mdpOptionProvider, "Applied-forces modules all implement MDP options");
        mdpOptionProvider->buildMdpOutput(builder);
    }

    IMDModule* fmmModule = impl_->modules_.at(std::string(FmmModuleInfo::sc_name)).get();
    const IMdpOptionProvider* fmmMdpOptionProvider = fmmModule->mdpOptionProvider();
    GMX_RELEASE_ASSERT(fmmMdpOptionProvider,
                       "Fast-multipole-method module must implement MDP options");
    fmmMdpOptionProvider->buildMdpOutput(builder);
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

ForceProviders* MDModules::initForceProviders(gmx_wallcycle* wallCycle)
{
    GMX_RELEASE_ASSERT(impl_->forceProviders_ == nullptr,
                       "Force providers initialized multiple times");
    impl_->forceProviders_ = std::make_unique<ForceProviders>(wallCycle);
    for (auto& [_, module] : impl_->modules_)
    {
        module->initForceProviders(impl_->forceProviders_.get());
    }
    return impl_->forceProviders_.get();
}

void MDModules::subscribeToPreProcessingNotifications()
{
    for (auto& [_, module] : impl_->modules_)
    {
        module->subscribeToPreProcessingNotifications(&impl_->notifiers_);
    }
}

void MDModules::subscribeToSimulationSetupNotifications()
{
    for (auto& [_, module] : impl_->modules_)
    {
        module->subscribeToSimulationSetupNotifications(&impl_->notifiers_);
    }
}

void MDModules::add(std::string_view nameView, std::shared_ptr<IMDModule> module)
{
    const std::string name(nameView);
    if (impl_->modules_.find(name) != impl_->modules_.end())
    {
        GMX_THROW(APIError(formatString("Module with name %s already added", name.c_str())));
    }
    impl_->modules_[name] = std::move(module);
}

const MDModulesNotifiers& MDModules::notifiers()
{
    return impl_->notifiers_;
}

} // namespace gmx
