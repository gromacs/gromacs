/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Implements Colvars MDModule class
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "colvarsMDModule.h"

#include <functional>
#include <memory>
#include <string>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "colvarsforceprovider.h"
#include "colvarsoptions.h"
#include "colvarssimulationsparameters.h"

enum class PbcType : int;
struct gmx_mtop_t;


namespace gmx
{
class IMDOutputProvider;
class IMdpOptionProvider;
class KeyValueTreeObject;
class MDLogger;

namespace
{

/*! \internal
 * \brief Colvars module
 *
 * Class that implements the colvars MDModule
 */
class ColvarsMDModule final : public IMDModule
{
public:
    //! \brief Construct the colvars module.
    explicit ColvarsMDModule() = default;

    // Now callbacks for several kinds of MdModuleNotification are created
    // and subscribed, and will be dispatched correctly at run time
    // based on the type of the parameter required by the lambda.

    /*! \brief Requests to be notified during pre-processing.
     *
     * \param[in] notifier allows the module to subscribe to notifications from MdModules.
     *
     * The Colvars MDModule subscribes to these notifications:
     *   - storing its internal parameters in a tpr file by writing to a
     *     key-value-tree during pre-processing by a function taking a
     *     KeyValueTreeObjectBuilder as parameter
     *   - Acess topology using gmx_mtop_t notification
     *   - Access MDLogger for notifications output
     *   - Access warninp for for grompp warnings output
     *   - Coordinates, PBC and box for setting up the proxy
     */
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifier) override
    {

        if (!colvarsOptions_.isActive())
        {
            return;
        }

        // Writing internal parameters during pre-processing
        const auto writeInternalParametersFunction = [this](KeyValueTreeObjectBuilder treeBuilder) {
            colvarsOptions_.writeInternalParametersToKvt(treeBuilder);
        };
        notifier->preProcessingNotifier_.subscribe(writeInternalParametersFunction);

        // Access of the topology during pre-processing
        const auto processTopologyFunction = [this](gmx_mtop_t* mtop) {
            colvarsOptions_.processTopology(mtop);
        };
        notifier->preProcessingNotifier_.subscribe(processTopologyFunction);

        // Set Logger during pre-processing
        const auto setLoggerFunction = [this](const MDLogger& logger) {
            colvarsOptions_.setLogger(logger);
        };
        notifier->preProcessingNotifier_.subscribe(setLoggerFunction);

        //  Notification of the Coordinates, box and pbc during pre-processing
        const auto processCoordinatesFunction = [this](const CoordinatesAndBoxPreprocessed& coord) {
            colvarsOptions_.processCoordinates(coord);
        };
        notifier->preProcessingNotifier_.subscribe(processCoordinatesFunction);

        //  Notification for the temperature
        const auto processTemperatureFunction = [this](const EnsembleTemperature& temp) {
            colvarsOptions_.processTemperature(temp);
        };
        notifier->preProcessingNotifier_.subscribe(processTemperatureFunction);
    }


    /*! \brief Request to be notified.
     * The Colvars MDModule subscribes to these notifications:
     *   - the LocalAtomSetManager sets in the simulation parameter setup
     *     by taking a LocalAtomSetManager * as parameter
     *   - the type of periodic boundary conditions that are used
     *     by taking a PeriodicBoundaryConditionType as parameter
     *   - the topology of the system
     *     by taking a gmx_mtop_t * as parameter
     *   - the communicator
     *     by taking a t_commrec as parameter
     *   - the simulation time step
     *     by taking a SimulationTimeStep as a parameter
     *   - MDLogger for notifications output
     *   - the .edr filename
     *     to provide a user-defined output prefix for Colvars output files
     */
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifier) override
    {
        if (!colvarsOptions_.isActive())
        {
            return;
        }

        // Reading internal parameters during simulation setup
        const auto readInternalParametersFunction = [this](const KeyValueTreeObject& tree) {
            colvarsOptions_.readInternalParametersFromKvt(tree);
        };
        notifier->simulationSetupNotifier_.subscribe(readInternalParametersFunction);
        // Retrieve the LocalAtomSetManager during simulation setup
        const auto setLocalAtomManagerFunction = [this](LocalAtomSetManager* localAtomSetManager) {
            this->ColvarsSimulationsParameters_.setLocalAtomSetManager(localAtomSetManager);
        };
        notifier->simulationSetupNotifier_.subscribe(setLocalAtomManagerFunction);

        // constructing PBC during simulation setup
        const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc) {
            this->ColvarsSimulationsParameters_.setPeriodicBoundaryConditionType(pbc);
        };
        notifier->simulationSetupNotifier_.subscribe(setPeriodicBoundaryContionsFunction);

        // Retrieve the topology during simulation setup
        const auto setTopologyFunction = [this](const gmx_mtop_t& mtop) {
            this->ColvarsSimulationsParameters_.setTopology(mtop);
        };
        notifier->simulationSetupNotifier_.subscribe(setTopologyFunction);

        // Retrieve the Communication Record during simulations setup
        const auto setCommFunction = [this](const t_commrec& cr) {
            this->ColvarsSimulationsParameters_.setComm(cr);
        };
        notifier->simulationSetupNotifier_.subscribe(setCommFunction);

        // setting the simulation time step
        const auto setSimulationTimeStepFunction = [this](const SimulationTimeStep& simulationTimeStep) {
            this->ColvarsSimulationsParameters_.setSimulationTimeStep(simulationTimeStep.delta_t);
        };
        notifier->simulationSetupNotifier_.subscribe(setSimulationTimeStepFunction);
        // Saving MDLogger during simulation setup
        const auto setLoggerFunction = [this](const MDLogger& logger) {
            this->ColvarsSimulationsParameters_.setLogger(logger);
        };
        notifier->simulationSetupNotifier_.subscribe(setLoggerFunction);

        const auto setEdrFileNameFunction = [this](const EdrOutputFilename& filename) {
            colvarsOptions_.processEdrFilename(filename);
        };
        notifier->simulationSetupNotifier_.subscribe(setEdrFileNameFunction);

        // writing checkpoint data
        const auto checkpointDataWriting = [this](MDModulesWriteCheckpointData checkpointData) {
            colvarsForceProvider_->writeCheckpointData(checkpointData, ColvarsModuleInfo::name_);
        };
        notifier->checkpointingNotifier_.subscribe(checkpointDataWriting);

        // reading checkpoint data
        const auto checkpointDataReading = [this](MDModulesCheckpointReadingDataOnMain checkpointData) {
            colvarsState_.readState(checkpointData.checkpointedData_, ColvarsModuleInfo::name_);
        };
        notifier->checkpointingNotifier_.subscribe(checkpointDataReading);

        // Handle the atoms redistributed signal
        const auto handleAtomsRedistributedSignal =
                [this](const MDModulesAtomsRedistributedSignal& atomsRedistributedSignal) {
                    colvarsForceProvider_->processAtomsRedistributedSignal(atomsRedistributedSignal);
                };
        notifier->simulationSetupNotifier_.subscribe(handleAtomsRedistributedSignal);
    }

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return &colvarsOptions_; }
    //! From IMDModule
    //! Colvars provide its own output
    IMDOutputProvider* outputProvider() override { return nullptr; }

    //! From IMDModule
    //! Add this module to the force providers if active
    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (colvarsOptions_.isActive())
        {

            colvarsForceProvider_ = std::make_unique<ColvarsForceProvider>(
                    colvarsOptions_.colvarsConfigContent(),
                    ColvarsSimulationsParameters_.topology(),
                    ColvarsSimulationsParameters_.periodicBoundaryConditionType(),
                    ColvarsSimulationsParameters_.logger(),
                    colvarsOptions_.colvarsInputFiles(),
                    colvarsOptions_.colvarsEnsTemp(),
                    colvarsOptions_.colvarsSeed(),
                    ColvarsSimulationsParameters_.localAtomSetManager(),
                    ColvarsSimulationsParameters_.comm(),
                    ColvarsSimulationsParameters_.simulationTimeStep(),
                    colvarsOptions_.colvarsAtomCoords(),
                    colvarsOptions_.colvarsOutputPrefix(),
                    colvarsState_);
            forceProviders->addForceProvider(colvarsForceProvider_.get());
        }
    }


private:
    //! The options provided for colvars
    ColvarsOptions colvarsOptions_;

    //! Parameters that become available at simulation setup time.
    ColvarsSimulationsParameters ColvarsSimulationsParameters_;

    //! Object that evaluates the forces
    std::unique_ptr<ColvarsForceProvider> colvarsForceProvider_;

    //! The state of colvars force provider to be written in the checkpoint
    ColvarsForceProviderState colvarsState_;
};


} // namespace

std::unique_ptr<IMDModule> ColvarsModuleInfo::create()
{
    return std::make_unique<ColvarsMDModule>();
}

const std::string ColvarsModuleInfo::name_ = "colvars";

} // namespace gmx
