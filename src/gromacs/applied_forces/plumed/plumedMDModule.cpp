/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Implements Plumed MDModule class
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "plumedMDModule.h"

#include <memory>
#include <string>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "plumedOptions.h"
#include "plumedforceprovider.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Plumed module
 *
 * Class that implements the plumed MDModule
 */
class PlumedMDModule final : public IMDModule
{
public:
    //! \brief Construct the plumed module.
    explicit PlumedMDModule() = default;
    // Now callbacks for several kinds of MdModuleNotification are created
    // and subscribed, and will be dispatched correctly at run time
    // based on the type of the parameter required by the lambda.

    /*! \brief Requests to be notified during pre-processing.
     *
     * Plumed does not act during the preprocessing phase of a simulation, so the input are ignored
     */
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /*notifier*/) override {}

    /*! \brief Subscribe to MDModules notifications for information needed just before the simulation.
     */
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifier) override
    {
        // TODO: add a check for threadmpi (see #5104, https://gitlab.com/gromacs/gromacs/-/merge_requests/4367#note_2102475958, the manual and the force provider for the details)

        // Access the plumed filename this is used to activate the plumed module
        notifier->simulationSetupNotifier_.subscribe(
                [this](const PlumedInputFilename& plumedFilename)
                { this->options_.setPlumedFile(plumedFilename.plumedFilename_); });
        // Access the temperature if it is constant during the simulation
        notifier->simulationSetupNotifier_.subscribe(
                [this](const EnsembleTemperature& ensembleT)
                { this->options_.setEnsembleTemperature(ensembleT); });
        // Access of the topology
        notifier->simulationSetupNotifier_.subscribe([this](const gmx_mtop_t& mtop)
                                                     { this->options_.setTopology(mtop); });
        // Retrieve the Communication Record during simulations setup
        notifier->simulationSetupNotifier_.subscribe([this](const t_commrec& cr)
                                                     { this->options_.setComm(cr); });
        // setting the simulation time step
        notifier->simulationSetupNotifier_.subscribe(
                [this](const SimulationTimeStep& simulationTimeStep)
                { this->options_.setSimulationTimeStep(simulationTimeStep.delta_t); });
        // Retrieve the starting behavior
        notifier->simulationSetupNotifier_.subscribe(
                [this](const StartingBehavior& startingBehavior)
                { this->options_.setStartingBehavior(startingBehavior); });
        //  writing checkpoint data
        notifier->checkpointingNotifier_.subscribe(
                [this](MDModulesWriteCheckpointData /*checkpointData*/)
                {
                    if (options_.active())
                    {
                        plumedForceProvider_->writeCheckpointData();
                    }
                });
    }

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    //! From IMDModule
    IMDOutputProvider* outputProvider() override
    { // Plumed provide its own output
        return nullptr;
    }
    //! From IMDModule - Adds this module to the force providers if active
    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (options_.active())
        {
            plumedForceProvider_ = std::make_unique<PlumedForceProvider>(options_.options());
            forceProviders->addForceProvider(plumedForceProvider_.get(), "Plumed");
        }
    }


private:
    //! Parameters that become available at simulation setup time.
    PlumedOptionProvider options_{};
    //! Object that evaluates the forces
    std::unique_ptr<PlumedForceProvider> plumedForceProvider_{};
};

} // namespace

std::unique_ptr<IMDModule> PlumedModuleInfo::create()
{
    return std::make_unique<PlumedMDModule>();
}
} // namespace gmx
