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

#include <string>

#include "gromacs/mdrun/isimulator.h"

#include "checkpointhelper.h"
#include "computeglobalselement.h"
#include "domdechelper.h"
#include "modularsimulatorinterfaces.h"
#include "pmeloadbalancehelper.h"

namespace gmx
{
class EnergyElement;
class EnergySignaller;
class FreeEnergyPerturbationElement;
class LoggingSignaller;
class ModularSimulator;
class NeighborSearchSignaller;
class ResetHandler;
class TopologyHolder;
class TrajectoryElementBuilder;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief The modular simulator
 *
 * Based on the input given, this simulator builds independent elements and
 * signallers and stores them in a respective vector. The run function
 * runs the simulation by, in turn, building a task list from the elements
 * for a predefined number of steps, then running the task list, and repeating
 * until the stop criterion is fulfilled.
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
    std::vector<SimulatorRunFunctionPtr> taskQueue_;
    //! The task iterator
    std::vector<SimulatorRunFunctionPtr>::const_iterator taskIterator_;

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
    std::vector<compat::not_null<ISimulatorElement*>> elementCallList_;

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
    SimulationSignals signals_;

    //! The topology
    std::unique_ptr<TopologyHolder> topologyHolder_;

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
        SignallerCallbackPtr registerLastStepCallback() override;
        //! INeighborSearchSignallerClient implementation
        SignallerCallbackPtr registerNSCallback() override;
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

/*!\internal
 * \brief Builder for ModularSimulatorAlgorithm objects
 *
 * TODO: The current builder automatically builds a simulator algorithm based on the
 *       input. This is only an intemediate step towards a builder that will create
 *       algorithms designed by the user of ModularSimulatorAlgorithm (for now, the
 *       only user is the ModularSimulator).
 * TODO: This mirrors all protected members of ISimulator. This hack allows to keep
 *       the number of line changes minimal, and will be removed as soon as the builder
 *       allows the user to compose the integrator algorithm.
 *       For the same reason, the constructElementsAndSignallers(), buildIntegrator(...),
 *       and buildForces(...) implementations were left in modularsimulator.cpp. See other
 *       to do - as the ModularSimulator will eventually design the algorithm, moving it
 *       would only cause unnecessary noise.
 */
class ModularSimulatorAlgorithmBuilder
{
public:
    //! Constructor
    explicit ModularSimulatorAlgorithmBuilder(compat::not_null<ModularSimulator*> simulator);
    //! Build algorithm
    ModularSimulatorAlgorithm build();

private:
    /*! \brief The initialisation
     *
     * This builds all signallers and elements, and is responsible to put
     * them in the correct order.
     */
    ModularSimulatorAlgorithm constructElementsAndSignallers();

    /*! \brief Build the integrator part of the simulator
     *
     * This includes the force calculation, state propagation, constraints,
     * global computation, and the points during the process at which valid
     * micro state / energy states are found. Currently, buildIntegrator
     * knows about NVE md and md-vv algorithms.
     */
    std::unique_ptr<ISimulatorElement>
    buildIntegrator(SignallerBuilder<NeighborSearchSignaller>* neighborSearchSignallerBuilder,
                    SignallerBuilder<EnergySignaller>*         energySignallerBuilder,
                    SignallerBuilder<LoggingSignaller>*        loggingSignallerBuilder,
                    SignallerBuilder<TrajectorySignaller>*     trajectorySignallerBuilder,
                    std::vector<ICheckpointHelperClient*>*     checkpointClients,
                    CheckBondedInteractionsCallbackPtr*        checkBondedInteractionsCallback,
                    compat::not_null<StatePropagatorData*>     statePropagatorDataPtr,
                    compat::not_null<EnergyElement*>           energyElementPtr,
                    FreeEnergyPerturbationElement*             freeEnergyPerturbationElementPtr,
                    bool                                       hasReadEkinState,
                    TopologyHolder*                            topologyHolder,
                    SimulationSignals*                         signals);

    //! Build the force element - can be normal forces or shell / flex constraints
    std::unique_ptr<ISimulatorElement>
    buildForces(SignallerBuilder<NeighborSearchSignaller>* neighborSearchSignallerBuilder,
                SignallerBuilder<EnergySignaller>*         energySignallerBuilder,
                StatePropagatorData*                       statePropagatorDataPtr,
                EnergyElement*                             energyElementPtr,
                FreeEnergyPerturbationElement*             freeEnergyPerturbationElement,
                TopologyHolder*                            topologyHolder);

    //! \cond
    //! Helper function to add elements or signallers to the call list via raw pointer
    template<typename T, typename U>
    static void addToCallList(U* element, std::vector<compat::not_null<T*>>& callList);
    //! Helper function to add elements or signallers to the call list via non-null raw pointer
    template<typename T, typename U>
    static void addToCallList(compat::not_null<U*> element, std::vector<compat::not_null<T*>>& callList);
    //! Helper function to add elements or signallers to the call list via smart pointer
    template<typename T, typename U>
    static void addToCallList(std::unique_ptr<U>& element, std::vector<compat::not_null<T*>>& callList);
    /*! \brief Helper function to add elements or signallers to the call list
     *         and move the ownership to the ownership list
     */
    template<typename T, typename U>
    static void addToCallListAndMove(std::unique_ptr<U>                 element,
                                     std::vector<compat::not_null<T*>>& callList,
                                     std::vector<std::unique_ptr<T>>&   elementList);
    //! \endcond

    //! Compute globals communication period
    const int nstglobalcomm_;

    // ISimulator data
    //! Handles logging.
    FILE* fplog;
    //! Handles communication.
    t_commrec* cr;
    //! Coordinates multi-simulations.
    const gmx_multisim_t* ms;
    //! Handles logging.
    const MDLogger& mdlog;
    //! Count of input file options.
    int nfile;
    //! Content of input file options.
    const t_filenm* fnm;
    //! Handles writing text output.
    const gmx_output_env_t* oenv;
    //! Contains command-line options to mdrun.
    const MdrunOptions& mdrunOptions;
    //! Whether the simulation will start afresh, or restart with/without appending.
    const StartingBehavior startingBehavior;
    //! Handles virtual sites.
    VirtualSitesHandler* vsite;
    //! Handles constraints.
    Constraints* constr;
    //! Handles enforced rotation.
    gmx_enfrot* enforcedRotation;
    //! Handles box deformation.
    BoxDeformation* deform;
    //! Handles writing output files.
    IMDOutputProvider* outputProvider;
    //! Handles notifications to MdModules for checkpoint writing
    const MdModulesNotifier& mdModulesNotifier;
    //! Contains user input mdp options.
    t_inputrec* inputrec;
    //! The Interactive Molecular Dynamics session.
    ImdSession* imdSession;
    //! The pull work object.
    pull_t* pull_work;
    //! The coordinate-swapping session.
    t_swap* swap;
    //! Full system topology.
    const gmx_mtop_t* top_global;
    //! Full simulation state (only non-nullptr on master rank).
    t_state* state_global;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory;
    //! Atom parameters for this domain.
    MDAtoms* mdAtoms;
    //! Manages flop accounting.
    t_nrnb* nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle;
    //! Parameters for force calculations.
    t_forcerec* fr;
    //! Data for energy output.
    gmx_enerdata_t* enerd;
    //! Kinetic energy data.
    gmx_ekindata_t* ekind;
    //! Schedule of work for each MD step for this task.
    MdrunScheduleWorkload* runScheduleWork;
    //! Parameters for replica exchange algorihtms.
    const ReplicaExchangeParameters& replExParams;
    //! Parameters for membrane embedding.
    gmx_membed_t* membed;
    //! Manages wall time accounting.
    gmx_walltime_accounting* walltime_accounting;
    //! Registers stop conditions
    StopHandlerBuilder* stopHandlerBuilder;
    //! Whether we're doing a rerun.
    bool doRerun;
};

//! \cond
template<typename T, typename U>
void ModularSimulatorAlgorithmBuilder::addToCallList(U* element, std::vector<compat::not_null<T*>>& callList)
{
    if (element)
    {
        callList.emplace_back(element);
    }
}

template<typename T, typename U>
void ModularSimulatorAlgorithmBuilder::addToCallList(gmx::compat::not_null<U*>          element,
                                                     std::vector<compat::not_null<T*>>& callList)
{
    callList.emplace_back(element);
}

template<typename T, typename U>
void ModularSimulatorAlgorithmBuilder::addToCallList(std::unique_ptr<U>&                element,
                                                     std::vector<compat::not_null<T*>>& callList)
{
    if (element)
    {
        callList.emplace_back(compat::make_not_null(element.get()));
    }
}

template<typename T, typename U>
void ModularSimulatorAlgorithmBuilder::addToCallListAndMove(std::unique_ptr<U> element,
                                                            std::vector<compat::not_null<T*>>& callList,
                                                            std::vector<std::unique_ptr<T>>& elementList)
{
    if (element)
    {
        callList.emplace_back(compat::make_not_null(element.get()));
        elementList.emplace_back(std::move(element));
    }
}
//! \endcond

} // namespace gmx

#endif // GROMACS_MODULARSIMULATOR_SIMULATORALGORITHM_H
