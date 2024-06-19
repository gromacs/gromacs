/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include <cstdio>

#include <any>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdrun/isimulator.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "checkpointhelper.h"
#include "domdechelper.h"
#include "freeenergyperturbationdata.h"
#include "modularsimulatorinterfaces.h"
#include "pmeloadbalancehelper.h"
#include "signallers.h"
#include "topologyholder.h"
#include "trajectoryelement.h"

struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_nrnb;

namespace gmx
{
enum class IntegrationStage;
class EnergyData;
class ModularSimulator;
class ResetHandler;
template<IntegrationStage integrationStage>
class Propagator;
class TopologyHolder;
class MDLogger;
class ReadCheckpointDataHolder;
class StatePropagatorData;
struct MdrunOptions;

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
                              const t_inputrec*        inputrec,
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
    //! List of schedulerElements (run calling sequence)
    std::vector<ISimulatorElement*> elementCallList_;
    //! List of schedulerElements (setup / teardown calling sequence)
    std::vector<ISimulatorElement*> elementSetupTeardownList_;
    //! List of pre-step scheduling functions
    std::vector<SchedulingFunction> preStepScheduling_;
    //! List of post-step scheduling functions
    std::vector<SchedulingFunction> postStepScheduling_;

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
    //! Arbitrary data with lifetime equal to the simulation (used to share data between elements)
    std::map<std::string, std::unique_ptr<std::any>> simulationData_;

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
    FILE* fpLog_;
    //! Handles communication.
    t_commrec* cr_;
    //! Handles logging.
    const MDLogger& mdLog_;
    //! Contains command-line options to mdrun.
    const MdrunOptions& mdrunOptions_;
    //! Contains user input mdp options.
    const t_inputrec* inputRec_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wallCycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
    //! Manages wall time accounting.
    gmx_walltime_accounting* wallTimeAccounting_;
};

/*! \internal
 * \brief Helper container with data connected to global communication
 *
 * This includes data that needs to be shared between elements involved in
 * global communication. This will become obsolete as soon as global
 * communication is moved to a client system (#3421 and #3887).
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

private:
    //! Compute globals communication period
    const int nstglobalcomm_;
    //! Signal vector (used by stop / reset / checkpointing signaller)
    SimulationSignals* simulationSignals_;
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
    template<typename Element>
    Element* storeElement(std::unique_ptr<Element> element);
    //! Check if an element is stored in the ModularSimulatorAlgorithmBuilder
    bool elementIsStored(const ISimulatorElement* element) const;
    /*! \brief Register callback to schedule a pre-step run
     *
     * This allows elements to schedule a function call before the integration step.
     * The function call is guaranteed to happen before any functions scheduled for
     * the integration step. It is not guaranteed to happen in any specific order
     * compared to other elements registering a pre-step scheduling function.
     */
    void registerPreStepScheduling(SchedulingFunction schedulingFunction);
    /*! \brief Register callback to schedule a post-step run
     *
     * This allows elements to schedule a function call after the integration step.
     * The function call is guaranteed to happen after all functions scheduled for
     * the integration step. It is not guaranteed to happen in any specific order
     * compared to other elements registering a post-step scheduling function.
     */
    void registerPostStepScheduling(SchedulingFunction schedulingFunction);
    /*! \brief Set arbitrary data in the ModularSimulatorAlgorithmBuilder
     *
     * Allows to store arbitrary data with lifetime equal to the builder. Functionality is used
     * by stateful static builder functions.
     */
    template<typename ValueType>
    void storeBuilderData(const std::string& key, const ValueType& value);
    //! Get previously stored builder data. Returns std::nullopt if key is not found.
    std::optional<std::any> builderData(const std::string& key) const;
    //! \copydoc ModularSimulatorAlgorithmBuilder::storeSimulationData()
    template<typename ValueType>
    void storeSimulationData(const std::string& key, ValueType&& value);
    //! \copydoc ModularSimulatorAlgorithmBuilder::simulationData()
    template<typename ValueType>
    std::optional<ValueType*> simulationData(const std::string& key);
    //! Register temperature / pressure control algorithm to be matched with a propagator
    void registerTemperaturePressureControl(std::function<void(const PropagatorConnection&)> registrationFunction);
    //! Register a propagator to be used with a temperature / pressure control algorithm
    void registerPropagator(PropagatorConnection connectionData);
    //! Register for callback after an update to the reference temperature
    void registerReferenceTemperatureUpdate(ReferenceTemperatureCallback referenceTemperatureCallback);
    //! Get a callback to change reference temperature
    ReferenceTemperatureCallback changeReferenceTemperatureCallback();

private:
    //! Pointer to the associated ModularSimulatorAlgorithmBuilder
    ModularSimulatorAlgorithmBuilder* builder_;
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
    //! Arbitrary data with lifetime equal to the builder (used by stateful static builder functions)
    std::map<std::string, std::any> builderData_;
    //! Arbitrary data with lifetime equal to the simulation (used to share data between elements)
    std::map<std::string, std::unique_ptr<std::any>> simulationData_;

    //! Pointer to the LegacySimulatorData object
    compat::not_null<LegacySimulatorData*> legacySimulatorData_;

    // Helper objects
    //! Signal vector (used by stop / reset / checkpointing signaller)
    std::unique_ptr<SimulationSignals> signals_;
    //! Helper object passed to element factory functions
    ModularSimulatorAlgorithmBuilderHelper elementAdditionHelper_;
    //! Container for minor aspects of global computation data
    GlobalCommunicationHelper globalCommunicationHelper_;
    //! Coordinates reduction for observables
    ObservablesReducer observablesReducer_;

    /*! \brief Set arbitrary data in the ModularSimulatorAlgorithm
     *
     * Allows to store arbitrary data with lifetime equal to the simulator algorithm.
     * Functionality allows elements to share arbitrary data.
     */
    template<typename ValueType>
    void storeSimulationData(const std::string& key, ValueType&& value);

    /*! \brief Get previously stored simulation data.
     *
     * Returns std::nullopt if key is not found.
     */
    template<typename ValueType>
    std::optional<ValueType*> simulationData(const std::string& key);

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
     * This will also add the element to the setup / teardown list, and register
     * it with all applicable signallers and infrastructure objects.
     * Note that simply adding an element using this function will not call it
     * during the simulation - it needs to be added to the call list separately.
     * Also note that generally, users will want to add elements to the call list,
     * but it might not be practical to do this in the same order.
     *
     * \tparam Element  Type of the Element
     * \param element  A unique pointer to the element
     * \return  A non-owning (raw) pointer to the element for further usage
     */
    template<typename Element>
    Element* addElementToSimulatorAlgorithm(std::unique_ptr<Element> element);

    /*! \brief Register existing element to infrastructure
     *
     * This function adds existing elements to the setup / teardown list, and
     * registers them with all applicable signallers and infrastructure objects.
     * This is only permissible for elements owned directly by the builder or
     * indirectly through data objects. Before registering the element, the function
     * checks that the element is owned by the builder or a known object.
     *
     * \tparam Element  Type of the Element
     * \param element   A non-owning (raw) pointer to the element
     */
    template<typename Element>
    void registerExistingElement(Element* element);

    /*! \brief Check if element is owned by *this
     *
     * \param element  Pointer to the element
     * \return  Bool indicating whether element is owned by *this
     */
    [[nodiscard]] bool elementExists(const ISimulatorElement* element) const;

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
    //! List of pre-step scheduling functions
    std::vector<SchedulingFunction> preStepScheduling_;
    //! List of post-step scheduling functions
    std::vector<SchedulingFunction> postStepScheduling_;

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
    //! Builder for the DomDecHelper
    DomDecHelperBuilder domDecHelperBuilder_;

    /*! \brief List of clients for the CheckpointHelper
     *
     * \todo Replace this by proper builder (#3422)
     */
    std::vector<ICheckpointHelperClient*> checkpointClients_;

    //! List of data to connect propagators to thermostats / barostats
    std::vector<PropagatorConnection> propagatorConnections_;
    //! List of temperature / pressure control registration functions
    std::vector<std::function<void(const PropagatorConnection&)>> pressureTemperatureControlRegistrationFunctions_;
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
 *             GlobalCommunicationHelper*              globalCommunicationHelper,
 *             ObservablesReducer*                     observablesReducer)
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
 * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
 * \param observablesReducer          Pointer to the \c ObservablesReducer object
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
                                     ObservablesReducer*         observablesReducer,
                                     Args&&... args)
{
    return Element::getElementPointerImpl(legacySimulatorData,
                                          builderHelper,
                                          statePropagatorData,
                                          energyData,
                                          freeEnergyPerturbationData,
                                          globalCommunicationHelper,
                                          observablesReducer,
                                          std::forward<Args>(args)...);
}

template<typename Element, typename... Args>
void ModularSimulatorAlgorithmBuilder::add(Args&&... args)
{
    if (algorithmHasBeenBuilt_)
    {
        GMX_THROW(SimulationAlgorithmSetupError(
                "Tried to add an element after ModularSimulationAlgorithm was built."));
    }

    // Get element from factory method
    auto* element = static_cast<Element*>(getElementPointer<Element>(legacySimulatorData_,
                                                                     &elementAdditionHelper_,
                                                                     statePropagatorData_.get(),
                                                                     energyData_.get(),
                                                                     freeEnergyPerturbationData_.get(),
                                                                     &globalCommunicationHelper_,
                                                                     &observablesReducer_,
                                                                     std::forward<Args>(args)...));

    // Make sure returned element pointer is owned by *this
    // Ensuring this makes sure we can control the life time
    if (!elementExists(element))
    {
        GMX_THROW(ElementNotFoundError("Tried to append non-existing element to call list."));
    }
    // Add to call list
    callList_.emplace_back(element);
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
    // Register element to DomDecHelper builder (if applicable)
    domDecHelperBuilder_.registerClient(castOrNull<IDomDecHelperClient, Element>(element));
}

template<typename Element>
Element* ModularSimulatorAlgorithmBuilder::addElementToSimulatorAlgorithm(std::unique_ptr<Element> element)
{
    // Store element
    elements_.emplace_back(std::move(element));
    // Get non-owning pointer for further use
    Element* elementPtr = static_cast<Element*>(elements_.back().get());
    // Register element to infrastructure
    registerExistingElement(elementPtr);

    return elementPtr;
}

template<typename Element>
void ModularSimulatorAlgorithmBuilder::registerExistingElement(Element* element)
{
    // Make sure the element pointer is owned by *this
    // Ensuring this makes sure we can control the life time
    if (!elementExists(element))
    {
        GMX_THROW(
                ElementNotFoundError("Tried to register non-existing element to infrastructure."));
    }

    // Add to setup / teardown list
    setupAndTeardownList_.emplace_back(element);
    // Register element to all applicable signallers
    registerWithInfrastructureAndSignallers(element);
}

template<typename Element>
Element* ModularSimulatorAlgorithmBuilderHelper::storeElement(std::unique_ptr<Element> element)
{
    return builder_->addElementToSimulatorAlgorithm(std::move(element));
}

template<typename ValueType>
void ModularSimulatorAlgorithmBuilderHelper::storeBuilderData(const std::string& key, const ValueType& value)
{
    builder_->builderData_[key] = std::any(value);
}

template<typename ValueType>
void ModularSimulatorAlgorithmBuilderHelper::storeSimulationData(const std::string& key, ValueType&& value)
{
    builder_->storeSimulationData(key, std::forward<ValueType>(value));
}

template<typename ValueType>
std::optional<ValueType*> ModularSimulatorAlgorithmBuilderHelper::simulationData(const std::string& key)
{
    return builder_->simulationData<ValueType>(key);
}

template<typename ValueType>
void ModularSimulatorAlgorithmBuilder::storeSimulationData(const std::string& key, ValueType&& value)
{
    GMX_RELEASE_ASSERT(simulationData_.count(key) == 0,
                       formatString("Key %s was already stored in simulation data.", key.c_str()).c_str());
    simulationData_[key] = std::make_unique<std::any>(std::forward<ValueType>(value));
    auto* ptrToData      = simulationData<ValueType>(key).value();
    registerWithInfrastructureAndSignallers(ptrToData);
}

template<typename ValueType>
std::optional<ValueType*> ModularSimulatorAlgorithmBuilder::simulationData(const std::string& key)
{
    const auto iter = simulationData_.find(key);
    if (iter == simulationData_.end())
    {
        return std::nullopt;
    }
    ValueType* data = std::any_cast<ValueType>(iter->second.get());
    GMX_RELEASE_ASSERT(data != nullptr,
                       formatString("Object stored in simulation data under key %s does not have "
                                    "the expected type.",
                                    key.c_str())
                               .c_str());
    return data;
}


} // namespace gmx

#endif // GROMACS_MODULARSIMULATOR_SIMULATORALGORITHM_H
