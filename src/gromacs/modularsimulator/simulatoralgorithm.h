/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Provides the modular simulator algorithm.
 *
 * Defines the ModularSimulatorAlgorithm class and its builder.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module.
 * Moving forward, the ModularSimulatorAlgorithmBuilder could be exposed to allow users to
 * create custom algorithm - currently algorithms are only created an used by the ModularSimulator,
 * meaning that this file is not exposed outside of the modular simulator module.
 */
#ifndef GROMACS_MODULARSIMULATOR_SIMULATORALGORITHM_H
#define GROMACS_MODULARSIMULATOR_SIMULATORALGORITHM_H

#include <any>
#include <map>
#include <optional>
#include <string>
#include <typeinfo>

#include "gromacs/mdrun/isimulator.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/exceptions.h"

#include "checkpointhelper.h"
#include "domdechelper.h"
#include "freeenergyperturbationdata.h"
#include "modularsimulatorinterfaces.h"
#include "pmeloadbalancehelper.h"
#include "signallers.h"
#include "topologyholder.h"
#include "trajectoryelement.h"

namespace gmx
{
enum class IntegrationStep;
class EnergyData;
class ModularSimulator;
class ResetHandler;
template<IntegrationStep algorithm>
class Propagator;
class TopologyHolder;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief The modular simulator
 *
 * Based on the input given, this simulator builds independent elements and
 * signallers and stores them in a respective vector. The run function
 * runs the simulation by, in turn, building a task list from the elements
 * for a predefined number of steps, then running the task list, and repeating
 * until the stop criterion is fulfilled.
 *
 * The simulator algorithm is owning all elements involved in the simulation
 * and is hence controlling their lifetime. This ensures that pointers and
 * callbacks exchanged between elements remain valid throughout the duration
 * of the simulation run.
 */
class ModularSimulatorAlgorithm final
{
public:
    //! Get next task in queue
    [[nodiscard]] const SimulatorRunFunction* getNextTask();

    // Only builder can construct
    friend class ModularSimulatorAlgorithmBuilder;

private:
    //! Constructor
    ModularSimulatorAlgorithm(std::string              topologyName,
                              FILE*                    fplog,
                              t_commrec*               cr,
                              const MDLogger&          mdlog,
                              const MdrunOptions&      mdrunOptions,
                              t_inputrec*              inputrec,
                              t_nrnb*                  nrnb,
                              gmx_wallcycle*           wcycle,
                              t_forcerec*              fr,
                              gmx_walltime_accounting* walltime_accounting);

    //! Calls setup for the simulator and all elements and signallers
    void setup();
    //! Calls teardown for the simulator and all elements
    void teardown();
    //! Re-populate task queue
    void updateTaskQueue();

    /*! \brief A function called once before the simulation run
     *
     * This function collects calls which need to be made once at the
     * beginning of the simulation run. These are mainly printing information
     * to the user and starting the wall time.
     */
    void simulatorSetup();

    /*! \brief A function called once after the simulation run
     *
     * This function collects calls which need to be made once at the
     * end of the simulation run.
     */
    void simulatorTeardown();

    /*! \brief A function called before every step
     *
     * This function collects calls which need to be made before every
     * simulation step. The need for this function could likely be
     * eliminated by rewriting the stop and reset handler to fit the
     * modular simulator approach.
     */
    void preStep(Step step, Time time, bool isNeighborSearchingStep);

    /*! \brief A function called after every step
     *
     * This function collects calls which need to be made after every
     * simulation step.
     */
    void postStep(Step step, Time time);

    /*! \brief Add run functions to the task queue
     *
     * This function calls PME load balancing and domain decomposition first,
     * and then queries the elements for all steps of the life time of the
     * current neighbor list for their run functions. When the next step
     * would be a neighbor-searching step, this function returns, such that
     * the simulator can run the pre-computed task list before updating the
     * neighbor list and re-filling the task list.
     *
     * TODO: In the current approach, unique pointers to tasks are created
     *       repeatedly. While preliminary measurements do not seem to indicate
     *       performance problems, a pool allocator would be an obvious
     *       optimization possibility, as would reusing the task queue if
     *       the task queue is periodic around the rebuild point (i.e. the
     *       task queue is identical between rebuilds).
     */
    void populateTaskQueue();

    //! The run queue
    std::vector<SimulatorRunFunction> taskQueue_;
    //! The task iterator
    std::vector<SimulatorRunFunction>::const_iterator taskIterator_;

    /* Note that the Simulator is owning the signallers and elements.
     * The ownership list and the call list of the elements are kept
     * separate, to allow to have elements more than once in the call
     * lists - for example, using velocity verlet, the compute globals
     * element needs to be scheduled more than once per step. For the
     * signallers, no distinction between ownership and call list is
     * made, all signallers are called exactly once per scheduling step.
     *
     * Objects being both a signaller and an element are not supported.
     */
    //! List of signalers (ownership and calling sequence)
    std::vector<std::unique_ptr<ISignaller>> signallerList_;
    //! List of schedulerElements (ownership)
    std::vector<std::unique_ptr<ISimulatorElement>> elementsOwnershipList_;
    //! List of schedulerElements (calling sequence)
    std::vector<ISimulatorElement*> elementCallList_;

    // Infrastructure elements
    //! The domain decomposition element
    std::unique_ptr<DomDecHelper> domDecHelper_;
    //! The PME load balancing element
    std::unique_ptr<PmeLoadBalanceHelper> pmeLoadBalanceHelper_;
    //! The checkpoint helper
    std::unique_ptr<CheckpointHelper> checkpointHelper_;
    //! The stop handler
    std::unique_ptr<StopHandler> stopHandler_;
    //! The reset handler
    std::unique_ptr<ResetHandler> resetHandler_;
    //! Signal vector (used by stop / reset / checkpointing signaller)
    std::unique_ptr<SimulationSignals> signals_;
    //! The topology
    std::unique_ptr<TopologyHolder> topologyHolder_;

    // Data structures
    //! The state propagator data
    std::unique_ptr<StatePropagatorData> statePropagatorData_;
    //! The energy data
    std::unique_ptr<EnergyData> energyData_;
    //! The free energy data
    std::unique_ptr<FreeEnergyPerturbationData> freeEnergyPerturbationData_;

    //! The current step
    Step step_;
    //! Whether the simulation run is finished
    bool runFinished_;

    /*! \internal
     * \brief Signal helper
     *
     * The simulator needs to know about the last step and the
     * neighbor searching step, which are determined in signallers.
     * This object can be registered to the signals and accessed by
     * the methods of the simulator.
     */
    class SignalHelper : public ILastStepSignallerClient, public INeighborSearchSignallerClient
    {
    public:
        //! The last step
        Step lastStep_ = std::numeric_limits<Step>::max();
        //! The next NS step
        Step nextNSStep_ = -1;
        //! ILastStepSignallerClient implementation
        std::optional<SignallerCallback> registerLastStepCallback() override;
        //! INeighborSearchSignallerClient implementation
        std::optional<SignallerCallback> registerNSCallback() override;
    };
    //! The signal helper object
    std::unique_ptr<SignalHelper> signalHelper_;

    // TODO: This is a hack for stop handler - needs to go once StopHandler
    //       is adapted to the modular simulator
    //! Whether this is a neighbor-searching step
    bool stophandlerIsNSStep_ = false;
    //! The current step
    Step stophandlerCurrentStep_ = -1;

    // ISimulator data - used for setup, teardown, pre step and post step
    //! Name of the topology
    std::string topologyName_;
    //! Handles logging.
    FILE* fplog;
    //! Handles communication.
    t_commrec* cr;
    //! Handles logging.
    const MDLogger& mdlog;
    //! Contains command-line options to mdrun.
    const MdrunOptions& mdrunOptions;
    //! Contains user input mdp options.
    t_inputrec* inputrec;
    //! Manages flop accounting.
    t_nrnb* nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle;
    //! Parameters for force calculations.
    t_forcerec* fr;
    //! Manages wall time accounting.
    gmx_walltime_accounting* walltime_accounting;
};

/*! \internal
 * \brief Helper container with data connected to global communication
 *
 * This includes data that needs to be shared between elements involved in
 * global communication. This will become obsolete as soon as global
 * communication is moved to a client system (#3421).
 */
class GlobalCommunicationHelper
{
public:
    //! Constructor
    GlobalCommunicationHelper(int nstglobalcomm, SimulationSignals* simulationSignals);

    //! Get the compute globals communication period
    [[nodiscard]] int nstglobalcomm() const;
    //! Get a pointer to the signals vector
    [[nodiscard]] SimulationSignals* simulationSignals();

    //! Set the callback to check the number of bonded interactions
    void setCheckBondedInteractionsCallback(CheckBondedInteractionsCallback callback);
    //! Move the callback to check the number of bonded interactions
    [[nodiscard]] CheckBondedInteractionsCallback moveCheckBondedInteractionsCallback();

private:
    //! Compute globals communication period
    const int nstglobalcomm_;
    //! Signal vector (used by stop / reset / checkpointing signaller)
    SimulationSignals* simulationSignals_;
    //! Callback to check the number of bonded interactions
    std::optional<CheckBondedInteractionsCallback> checkBondedInteractionsCallback_;
};

class ModularSimulatorAlgorithmBuilder;

/*! \internal
 * \brief Helper for element addition
 *
 * Such an object will be given to each invocation of getElementPointer
 *
 * Note: It would be nicer to define this as a member type of
 * ModularSimulatorAlgorithmBuilder, but this would break forward declaration.
 * This object is therefore defined as friend class.
 */
class ModularSimulatorAlgorithmBuilderHelper
{
public:
    //! Constructor
    ModularSimulatorAlgorithmBuilderHelper(ModularSimulatorAlgorithmBuilder* builder);
    //! Store an element to the ModularSimulatorAlgorithmBuilder
    ISimulatorElement* storeElement(std::unique_ptr<ISimulatorElement> element);
    //! Check if an element is stored in the ModularSimulatorAlgorithmBuilder
    bool elementIsStored(const ISimulatorElement* element) const;
    //! Set arbitrary data in the ModularSimulatorAlgorithmBuilder. Helpful for stateful elements.
    template<typename ValueType>
    void storeValue(const std::string& key, const ValueType& value);
    //! Get previously stored data. Returns std::nullopt if key is not found.
    std::optional<std::any> getStoredValue(const std::string& key) const;
    //! Register a thermostat that accepts propagator registrations
    void registerThermostat(std::function<void(const PropagatorThermostatConnection&)> registrationFunction);
    //! Register a barostat that accepts propagator registrations
    void registerBarostat(std::function<void(const PropagatorBarostatConnection&)> registrationFunction);
    //! Register a propagator to the thermostat used
    void registerWithThermostat(PropagatorThermostatConnection connectionData);
    //! Register a propagator to the barostat used
    void registerWithBarostat(PropagatorBarostatConnection connectionData);

private:
    //! Pointer to the associated ModularSimulatorAlgorithmBuilder
    ModularSimulatorAlgorithmBuilder* builder_;
    std::map<std::string, std::any>   values_;
};

/*!\internal
 * \brief Builder for ModularSimulatorAlgorithm objects
 *
 * This builds a ModularSimulatorAlgorithm.
 *
 * Users can add elements and define their call order by calling the templated
 * add<Element> function. Note that only elements that have a static
 * getElementPointerImpl factory method can be built in that way.
 *
 * Note that each ModularSimulatorAlgorithmBuilder can only be used to build
 * one ModularSimulatorAlgorithm object, i.e. build() can only be called once.
 * During the call to build, all elements and other infrastructure objects will
 * be moved to the built ModularSimulatorAlgorithm object, such that further use
 * of the builder would not make sense.
 * Any access to the build or add<> methods after the first call to
 * build() will result in an exception being thrown.
 */
class ModularSimulatorAlgorithmBuilder final
{
public:
    //! Constructor
    ModularSimulatorAlgorithmBuilder(compat::not_null<LegacySimulatorData*>    legacySimulatorData,
                                     std::unique_ptr<ReadCheckpointDataHolder> checkpointDataHolder);
    //! Build algorithm
    ModularSimulatorAlgorithm build();

    /*! \brief  Add element to the modular simulator algorithm builder
     *
     * This function has a general implementation, which will call the getElementPointer(...)
     * factory function.
     *
     * \tparam Element  The element type
     * \tparam Args     A variable number of argument types
     * \param args      A variable number of arguments
     */
    template<typename Element, typename... Args>
    void add(Args&&... args);

    //! Allow access from helper
    friend class ModularSimulatorAlgorithmBuilderHelper;

private:
    //! The state of the builder
    bool algorithmHasBeenBuilt_ = false;

    // Data structures
    //! The state propagator data
    std::unique_ptr<StatePropagatorData> statePropagatorData_;
    //! The energy data
    std::unique_ptr<EnergyData> energyData_;
    //! The free energy data
    std::unique_ptr<FreeEnergyPerturbationData> freeEnergyPerturbationData_;

    //! Pointer to the LegacySimulatorData object
    compat::not_null<LegacySimulatorData*> legacySimulatorData_;

    // Helper objects
    //! Signal vector (used by stop / reset / checkpointing signaller)
    std::unique_ptr<SimulationSignals> signals_;
    //! Helper object passed to element factory functions
    ModularSimulatorAlgorithmBuilderHelper elementAdditionHelper_;
    //! Container for global computation data
    GlobalCommunicationHelper globalCommunicationHelper_;

    /*! \brief  Register an element to all applicable signallers and infrastructure elements
     *
     * \tparam Element  Type of the Element
     * \param element   Pointer to the element
     */
    template<typename Element>
    void registerWithInfrastructureAndSignallers(Element* element);

    /*! \brief Take ownership of element
     *
     * This function returns a non-owning pointer to the new location of that
     * element, allowing further usage (e.g. adding the element to the call list).
     * Note that simply addin an element using this function will not call it
     * during the simulation - it needs to be added to the call list separately.
     * Note that generally, users will want to add elements to the call list, but
     * it might not be practical to do this in the same order.
     *
     * \param element  A unique pointer to the element
     * \return  A non-owning (raw) pointer to the element for further usage
     */
    ISimulatorElement* addElementToSimulatorAlgorithm(std::unique_ptr<ISimulatorElement> element);

    /*! \brief Check if element is owned by *this
     *
     * \param element  Pointer to the element
     * \return  Bool indicating whether element is owned by *this
     */
    [[nodiscard]] bool elementExists(const ISimulatorElement* element) const;

    /*! \brief Add element to setupAndTeardownList_ if it's not already there
     *
     * \param element  Element pointer to be added
     */
    void addElementToSetupTeardownList(ISimulatorElement* element);

    //! Vector to store elements, allowing the SimulatorAlgorithm to control their lifetime
    std::vector<std::unique_ptr<ISimulatorElement>> elements_;
    /*! \brief List defining in which order elements are called every step
     *
     * Elements may be referenced more than once if they should be called repeatedly
     */
    std::vector<ISimulatorElement*> callList_;
    /*! \brief  List defining in which order elements are set up and torn down
     *
     * Elements should only appear once in this list
     */
    std::vector<ISimulatorElement*> setupAndTeardownList_;

    //! Builder for the NeighborSearchSignaller
    SignallerBuilder<NeighborSearchSignaller> neighborSearchSignallerBuilder_;
    //! Builder for the LastStepSignaller
    SignallerBuilder<LastStepSignaller> lastStepSignallerBuilder_;
    //! Builder for the LoggingSignaller
    SignallerBuilder<LoggingSignaller> loggingSignallerBuilder_;
    //! Builder for the EnergySignaller
    SignallerBuilder<EnergySignaller> energySignallerBuilder_;
    //! Builder for the TrajectorySignaller
    SignallerBuilder<TrajectorySignaller> trajectorySignallerBuilder_;
    //! Builder for the TrajectoryElementBuilder
    TrajectoryElementBuilder trajectoryElementBuilder_;
    //! Builder for the TopologyHolder
    TopologyHolder::Builder topologyHolderBuilder_;
    //! Builder for the CheckpointHelper
    CheckpointHelperBuilder checkpointHelperBuilder_;

    /*! \brief List of clients for the CheckpointHelper
     *
     * \todo Replace this by proper builder (#3422)
     */
    std::vector<ICheckpointHelperClient*> checkpointClients_;

    //! List of thermostat registration functions
    std::vector<std::function<void(const PropagatorThermostatConnection&)>> thermostatRegistrationFunctions_;
    //! List of barostat registration functions
    std::vector<std::function<void(const PropagatorBarostatConnection&)>> barostatRegistrationFunctions_;
    //! List of data to connect propagators to thermostats
    std::vector<PropagatorThermostatConnection> propagatorThermostatConnections_;
    //! List of data to connect propagators to barostats
    std::vector<PropagatorBarostatConnection> propagatorBarostatConnections_;
};

/*! \internal
 * \brief Factory function for elements that can be added via ModularSimulatorAlgorithmBuilder:
 *        Get a pointer to an object of type \c Element to add to the call list
 *
 * This allows elements to be built via the templated ModularSimulatorAlgorithmBuilder::add<Element>
 * method. Elements buildable throught this factor function are required to implement a static
 * function with minimal signature
 *
 *     static ISimulatorElement* getElementPointerImpl(
 *             LegacySimulatorData*                    legacySimulatorData,
 *             ModularSimulatorAlgorithmBuilderHelper* builderHelper,
 *             StatePropagatorData*                    statePropagatorData,
 *             EnergyData*                             energyData,
 *             FreeEnergyPerturbationData*             freeEnergyPerturbationData,
 *             GlobalCommunicationHelper*              globalCommunicationHelper)
 *
 * This function may also accept additional parameters which are passed using the variadic
 * template parameter pack forwarded in getElementPointer.
 *
 * This function returns a pointer to an object of the Element type. Note that the caller will
 * check whether the returned object has previously been stored using the `storeElement`
 * function, and throw an exception if the element is not found.
 * The function can check whether a previously stored pointer is valid using
 * the `checkElementExistence` function. Most implementing functions will simply want
 * to create an object, store it using `storeElement`, and then use the return value of
 * `storeElement` as a return value to the caller. However, this setup allows the function
 * to store a created element (using a static pointer inside the function) and return it
 * in case that the factory function is called repeatedly. This allows to create an element
 * once, but have it called multiple times during the simulation run.
 *
 * \see ModularSimulatorAlgorithmBuilder::add
 *      Function using this functionality
 * \see ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>::getElementPointerImpl
 *      Implementation using the single object / multiple call sites functionality
 *
 * \tparam Element The type of the element
 * \tparam Args  Variable number of argument types allowing specific implementations to have
 *               additional arguments
 *
 * \param legacySimulatorData  Pointer allowing access to simulator level data
 * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
 * \param statePropagatorData  Pointer to the \c StatePropagatorData object
 * \param energyData  Pointer to the \c EnergyData object
 * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
 * \param globalCommunicationHelper  Pointer to the \c GlobalCommunicationHelper object
 * \param args  Variable number of additional parameters to be forwarded
 *
 * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
 */
template<typename Element, typename... Args>
ISimulatorElement* getElementPointer(LegacySimulatorData*                    legacySimulatorData,
                                     ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                     StatePropagatorData*                    statePropagatorData,
                                     EnergyData*                             energyData,
                                     FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                     GlobalCommunicationHelper*  globalCommunicationHelper,
                                     Args&&... args)
{
    return Element::getElementPointerImpl(legacySimulatorData, builderHelper, statePropagatorData,
                                          energyData, freeEnergyPerturbationData,
                                          globalCommunicationHelper, std::forward<Args>(args)...);
}

template<typename Element, typename... Args>
void ModularSimulatorAlgorithmBuilder::add(Args&&... args)
{
    if (algorithmHasBeenBuilt_)
    {
        throw SimulationAlgorithmSetupError(
                "Tried to add an element after ModularSimulationAlgorithm was built.");
    }

    // Get element from factory method
    auto* element = static_cast<Element*>(getElementPointer<Element>(
            legacySimulatorData_, &elementAdditionHelper_, statePropagatorData_.get(),
            energyData_.get(), freeEnergyPerturbationData_.get(), &globalCommunicationHelper_,
            std::forward<Args>(args)...));

    // Make sure returned element pointer is owned by *this
    // Ensuring this makes sure we can control the life time
    if (!elementExists(element))
    {
        throw ElementNotFoundError("Tried to append non-existing element to call list.");
    }
    // Add to call list
    callList_.emplace_back(element);
    // Add to setup / teardown list if element hasn't been added yet
    addElementToSetupTeardownList(element);
    // Register element to all applicable signallers
    registerWithInfrastructureAndSignallers(element);
}

//! Returns a pointer casted to type Base if the Element is derived from Base
template<typename Base, typename Element>
static std::enable_if_t<std::is_base_of<Base, Element>::value, Base*> castOrNull(Element* element)
{
    return static_cast<Base*>(element);
}

//! Returns a nullptr of type Base if Element is not derived from Base
template<typename Base, typename Element>
static std::enable_if_t<!std::is_base_of<Base, Element>::value, Base*> castOrNull(Element gmx_unused* element)
{
    return nullptr;
}

template<typename Element>
void ModularSimulatorAlgorithmBuilder::registerWithInfrastructureAndSignallers(Element* element)
{
    // Register element to all applicable signallers
    neighborSearchSignallerBuilder_.registerSignallerClient(
            castOrNull<INeighborSearchSignallerClient, Element>(element));
    lastStepSignallerBuilder_.registerSignallerClient(castOrNull<ILastStepSignallerClient, Element>(element));
    loggingSignallerBuilder_.registerSignallerClient(castOrNull<ILoggingSignallerClient, Element>(element));
    energySignallerBuilder_.registerSignallerClient(castOrNull<IEnergySignallerClient, Element>(element));
    trajectorySignallerBuilder_.registerSignallerClient(
            castOrNull<ITrajectorySignallerClient, Element>(element));
    // Register element to trajectory element (if applicable)
    trajectoryElementBuilder_.registerWriterClient(castOrNull<ITrajectoryWriterClient, Element>(element));
    // Register element to topology holder (if applicable)
    topologyHolderBuilder_.registerClient(castOrNull<ITopologyHolderClient, Element>(element));
    // Register element to checkpoint client (if applicable)
    checkpointHelperBuilder_.registerClient(castOrNull<ICheckpointHelperClient, Element>(element));
}


template<typename ValueType>
void ModularSimulatorAlgorithmBuilderHelper::storeValue(const std::string& key, const ValueType& value)
{
    values_[key] = std::any(value);
}


} // namespace gmx

#endif // GROMACS_MODULARSIMULATOR_SIMULATORALGORITHM_H
