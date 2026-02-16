/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Declares the FmmMDModule class, an implementation of the IMDModule interface.
 *
 * \author Muhammad Umair Sadiq <mumairsadiq1@gmail.com>
 */

#include "gmxpre.h"

#include "fmm_mdmodule.h"

#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"

#include "fmm_mdpoptions.h"
#include "fmm_mdpvalidator.h"
#include "fmmforceproviderbuilder.h"

namespace gmx
{

namespace
{
/*! \internal
 * \brief FMM module
 *
 * Class that implements the FMM MDModule
 */
class FmmMDModule final : public IMDModule
{
public:
    //! \brief Construct the FMM module.
    explicit FmmMDModule() = default;

    /*! \brief Subscribes to preprocessing (grompp) stage notifications.
     *
     * Registers handlers needed during preprocessing. This includes:
     *   - Setting the WarningHandler to emit preprocessing-time errors.
     *   - Checking that the Coulomb interaction type is set to FMM.
     *   - Validating FMM MDP settings.
     */
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override
    {
        if (fmmMdpOptions_.activeFmmBackend() == ActiveFmmBackend::Inactive)
        {
            return;
        }

        fmmMdpValidator_ = std::make_unique<FmmMdpValidator>(*fmmMdpOptions_.activeFmmOptions());

        const auto setWarningFunction = [this](WarningHandler* wi)
        { fmmMdpValidator_->setWarningHandler(wi); };
        notifiers->preProcessingNotifier_.subscribe(setWarningFunction);

        // Make sure that the Coulomb type is FMM and validate the FMM-related .mdp options
        const auto setCoulombTypeFunction = [this](const MdModulesCoulombTypeInfo& coulombTypeInfo)
        { fmmMdpValidator_->validateFmmMdpSettings(coulombTypeInfo); };
        notifiers->preProcessingNotifier_.subscribe(setCoulombTypeFunction);
    }

    /*! \brief Subscribes to simulation setup notifications.
     *
     * \param[in] notifiers Pointer to the simulation notifiers used to subscribe to setup events.
     *
     * The FMM module registers for the following simulation setup notifications:
     *   - gmx_mtop_t: Provides the topology, could be need for setting up exclusions in the FMM direct interactions.
     *   - PbcType : Supplies periodic boundary condition type used during simulation.
     *   - MDLogger: Provides access to the simulation logger for diagnostic or status messages.
     *   - MDModulesDirectProvider: Indicates if FMM handles short-range coulomb interactions.
     */
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override
    {
        fmmNotifiedDirectProvider_ = false;
        if (fmmMdpOptions_.activeFmmBackend() == ActiveFmmBackend::Inactive)
        {
            return;
        }

        fmmForceProviderBuilder_.subscribeToSimulationSetupNotifications(notifiers);


        FmmDirectProvider directProvider = fmmMdpOptions_.directProvider();
        // Subscribe to coulomb short range handler notification to indicate if FMM handles short-range coulomb interactions.
        const auto setCoulombDirectProviderFlag =
                [directProvider, this](MDModulesDirectProvider* mdModulesCoulombDirectProvider)
        {
            mdModulesCoulombDirectProvider->isDirectProvider = (directProvider == FmmDirectProvider::Fmm);
            fmmNotifiedDirectProvider_ = true;
        };
        notifiers->simulationSetupNotifier_.subscribe(setCoulombDirectProviderFlag);
    }

    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    void initForceProviders(ForceProviders* forceProviders) override
    {

        const ActiveFmmBackend fmmBackend = fmmMdpOptions_.activeFmmBackend();
        if (fmmBackend == ActiveFmmBackend::Inactive)
        {
            return;
        }

        if (!fmmNotifiedDirectProvider_)
        {
            GMX_THROW(MDModuleSetupError(
                    "FMM module must be notified about direct provider before force "
                    "provider initialization"));
        }

        fmmForceProviderBuilder_.setFmmOptions(fmmMdpOptions_.activeFmmOptions());
        fmmForceProvider_ = fmmForceProviderBuilder_.build();
        forceProviders->addForceProvider(fmmForceProvider_.get(), std::string(FmmModuleInfo::sc_name));
    }

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return &fmmMdpOptions_; }

    //! From IMDModule
    //! doesn't need extra output
    IMDOutputProvider* outputProvider() override { return nullptr; }

private:
    std::unique_ptr<FmmForceProvider> fmmForceProvider_;
    FmmMdpOptions                     fmmMdpOptions_;
    std::unique_ptr<FmmMdpValidator>  fmmMdpValidator_;
    FmmForceProviderBuilder           fmmForceProviderBuilder_;
    bool                              fmmNotifiedDirectProvider_ = false;
};

} // namespace

std::unique_ptr<IMDModule> FmmModuleInfo::create()
{
    return std::make_unique<FmmMDModule>();
}

} // namespace gmx
