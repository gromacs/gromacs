/*! \internal \file
 * \brief
 * Implements RAMD class that implements IMDModule interface
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "ramd.h"
#include "ramdoptions.h"
#include "ramdoutputprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief RAMD module
 */
class RAMD final : public IMDModule
{
public:
    //! \brief Construct the RAMD module.
    explicit RAMD() = default;

//     // Now callbacks for several kinds of MdModuleNotification are created
//     // and subscribed, and will be dispatched correctly at run time
//     // based on the type of the parameter required by the lambda.

//     /*! \brief Requests to be notified during pre-processing.
//      *
//      * \param[in] notifiers allows the module to subscribe to notifications from MdModules.
//      *
//      * The RAMD code subscribes to these notifications:
//      *   - setting atom group indices in the RAMDOptions_ from an
//      *     index group string by taking a parmeter const IndexGroupsAndNames &
//      *   - storing its internal parameters in a tpr file by writing to a
//      *     key-value-tree during pre-processing by a function taking a
//      *     KeyValueTreeObjectBuilder as parameter
//      *   - Modify topology according to RAMD rules using gmx_mtop_t notification
//      *     and utilizing RAMDTopologyPreprocessor class
//      *   - Access MDLogger for notifications output
//      *   - Access warninp for for grompp warnings output
//      *   - Coordinates, PBC and box for CP2K input generation
//      *   - QM Input file provided with -qmi option of grompp
//      */
//     void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override
//     {
//         if (!RAMDOptions_.active())
//         {
//             return;
//         }

//         // Writing internal parameters during pre-processing
//         const auto writeInternalParametersFunction = [this](KeyValueTreeObjectBuilder treeBuilder)
//         { RAMDOptions_.writeInternalParametersToKvt(treeBuilder); };
//         notifiers->preProcessingNotifier_.subscribe(writeInternalParametersFunction);

//         // Setting atom group indices
//         const auto setRAMDGroupIndicesFunction = [this](const IndexGroupsAndNames& indexGroupsAndNames)
//         { RAMDOptions_.setRAMDGroupIndices(indexGroupsAndNames); };
//         notifiers->preProcessingNotifier_.subscribe(setRAMDGroupIndicesFunction);

//         // Set Logger during pre-processing
//         const auto setLoggerFunction = [this](const MDLogger& logger)
//         { RAMDOptions_.setLogger(logger); };
//         notifiers->preProcessingNotifier_.subscribe(setLoggerFunction);

//         // Set warning output during pre-processing
//         const auto setWarninpFunction = [this](WarningHandler* wi) { RAMDOptions_.setWarninp(wi); };
//         notifiers->preProcessingNotifier_.subscribe(setWarninpFunction);

//         // Notification of the Coordinates, box and pbc during pre-processing
//         const auto processCoordinatesFunction = [this](const CoordinatesAndBoxPreprocessed& coord)
//         { RAMDOptions_.processCoordinates(coord); };
//         notifiers->preProcessingNotifier_.subscribe(processCoordinatesFunction);

//         // Modification of the topology during pre-processing
//         const auto modifyRAMDTopologyFunction = [this](gmx_mtop_t* mtop)
//         { RAMDOptions_.modifyRAMDTopology(mtop); };
//         notifiers->preProcessingNotifier_.subscribe(modifyRAMDTopologyFunction);

//         // Notification of the QM input file provided via -qmi option of grompp
//         const auto setQMExternalInputFileNameFunction = [this](const QMInputFileName& qmInputFileName)
//         { RAMDOptions_.setQMExternalInputFile(qmInputFileName); };
//         notifiers->preProcessingNotifier_.subscribe(setQMExternalInputFileNameFunction);
//     }

//     /*! \brief Requests to be notified during simulation setup.
//      * The RAMD code subscribes to these notifications:
//      *   - reading its internal parameters from a key-value-tree during
//      *     simulation setup by taking a const KeyValueTreeObject & parameter
//      *   - *.tpr filename for CP2K input generation
//      *   - constructing local atom sets in the simulation parameter setup
//      *     by taking a LocalAtomSetManager * as parameter
//      *   - the type of periodic boundary conditions that are used
//      *     by taking a PeriodicBoundaryConditionType as parameter
//      *   - Access MDLogger for notifications output
//      *   - Disable PME-only ranks for RAMD runs
//      *   - Request QM energy output to md.log
//      */
//     void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override
//     {
//         if (!RAMDOptions_.active())
//         {
//             return;
//         }

//         // Reading internal parameters during simulation setup
//         const auto readInternalParametersFunction = [this](const KeyValueTreeObject& tree)
//         { RAMDOptions_.readInternalParametersFromKvt(tree); };
//         notifiers->simulationSetupNotifier_.subscribe(readInternalParametersFunction);

//         // Process tpr filename
//         const auto setTprFileNameFunction = [this](const MdRunInputFilename& tprName)
//         { RAMDOptions_.processTprFilename(tprName); };
//         notifiers->simulationSetupNotifier_.subscribe(setTprFileNameFunction);

//         // constructing local atom sets during simulation setup
//         const auto setLocalAtomSetFunction = [this](LocalAtomSetManager* localAtomSetManager)
//         {
//             LocalAtomSet atomSet1 = localAtomSetManager->add(RAMDOptions_.parameters().qmIndices_);
//             this->RAMDSimulationParameters_.setLocalQMAtomSet(atomSet1);
//             LocalAtomSet atomSet2 = localAtomSetManager->add(RAMDOptions_.parameters().mmIndices_);
//             this->RAMDSimulationParameters_.setLocalMMAtomSet(atomSet2);
//         };
//         notifiers->simulationSetupNotifier_.subscribe(setLocalAtomSetFunction);

//         // Reading PBC parameters during simulation setup
//         const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc)
//         { this->RAMDSimulationParameters_.setPeriodicBoundaryConditionType(pbc); };
//         notifiers->simulationSetupNotifier_.subscribe(setPeriodicBoundaryContionsFunction);

//         // Saving MDLogger during simulation setup
//         const auto setLoggerFunction = [this](const MDLogger& logger)
//         { this->RAMDSimulationParameters_.setLogger(logger); };
//         notifiers->simulationSetupNotifier_.subscribe(setLoggerFunction);

//         // Adding output to energy file
//         const auto requestEnergyOutput = [](MDModulesEnergyOutputToRAMDRequestChecker* energyOutputRequest)
//         { energyOutputRequest->energyOutputToRAMD_ = true; };
//         notifiers->simulationSetupNotifier_.subscribe(requestEnergyOutput);

//         // Request to disable PME-only ranks, which are not compatible with CP2K
//         const auto requestPmeRanks = [](SeparatePmeRanksPermitted* pmeRanksPermitted)
//         {
//             pmeRanksPermitted->disablePmeRanks(
//                     "Separate PME-only ranks are not compatible with RAMD MdModule");
//         };
//         notifiers->simulationSetupNotifier_.subscribe(requestPmeRanks);
//     }

//     //! No subscriptions to run notifications
//     void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    //! Returns an interface for handling mdp input (and tpr I/O).
    IMdpOptionProvider* mdpOptionProvider() override { return &RAMDOptions_; }
    //! Returns an interface for handling output files during simulation.
    IMDOutputProvider* outputProvider() override { return &RAMDOutputProvider_; }
    //! Initializes force providers from this module.
    void initForceProviders(ForceProviders* forceProviders) override {}
    //! Subscribe to pre processing notifications
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override {}
    //! Subscribe to simulation setup notifications
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override {}
    //! Subscribe to simulation run notifications
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* notifiers) override {}

private:
    //! The output provider
    RAMDOutputProvider RAMDOutputProvider_;

    //! The options provided for RAMD
    RAMDOptions RAMDOptions_;

    //! Object that evaluates the forces
    // std::unique_ptr<RAMDForceProvider> forceProvider_;
    /*! \brief Parameters for RAMD that become available at
     * simulation setup time.
     */
    // RAMDSimulationParameterSetup RAMDSimulationParameters_;

    GMX_DISALLOW_COPY_AND_ASSIGN(RAMD);
};

} // namespace 

std::unique_ptr<IMDModule> RAMDModuleInfo::create()
{
    return std::make_unique<RAMD>();
}

} // namespace gmx
