/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Provides the modular simulator.
 *
 * Defines the ModularSimulator class. Provides checkUseModularSimulator() utility function
 * to determine whether the ModularSimulator should be used.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */
#ifndef GROMACS_MODULARSIMULATOR_MODULARSIMULATOR_H
#define GROMACS_MODULARSIMULATOR_MODULARSIMULATOR_H

#include <queue>

#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdrun/isimulator.h"

#include "checkpointhelper.h"
#include "domdechelper.h"
#include "modularsimulatorinterfaces.h"
#include "pmeloadbalancehelper.h"
#include "topologyholder.h"

namespace gmx
{
class DomDecHelper;
class EnergyElement;
class EnergySignaller;
class FreeEnergyPerturbationElement;
class LoggingSignaller;
class StatePropagatorData;
class NeighborSearchSignaller;
class PmeLoadBalanceHelper;
class TrajectoryElementBuilder;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief The modular simulator
 *
 * Based on the input given, this simulator builds independent elements and
 * signallers and stores them in a respective vector. The run function
 * runs the simulation by, in turn, building a task list from the elements
 * for a predefined number of steps, then running the task list, and repeating
 * until the stop criterion is fulfilled.
 */
class ModularSimulator final : public ISimulator
{
public:
    //! Run the simulator
    void run() override;

    //! Check for disabled functionality
    static bool isInputCompatible(bool                             exitOnFailure,
                                  const t_inputrec*                inputrec,
                                  bool                             doRerun,
                                  const gmx_mtop_t&                globalTopology,
                                  const gmx_multisim_t*            ms,
                                  const ReplicaExchangeParameters& replExParams,
                                  const t_fcdata*                  fcd,
                                  bool                             doEssentialDynamics,
                                  bool                             doMembed);

    // Only builder can construct
    friend class SimulatorBuilder;

private:
    //! Constructor
    template<typename... Args>
    explicit ModularSimulator(Args&&... args);

    /*! \brief The initialisation
     *
     * This builds all signallers and elements, and is responsible to put
     * them in the correct order.
     */
    void constructElementsAndSignallers();

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
                    TrajectoryElementBuilder*                  trajectoryElementBuilder,
                    std::vector<ICheckpointHelperClient*>*     checkpointClients,
                    CheckBondedInteractionsCallbackPtr*        checkBondedInteractionsCallback,
                    compat::not_null<StatePropagatorData*>     statePropagatorDataPtr,
                    compat::not_null<EnergyElement*>           energyElementPtr,
                    FreeEnergyPerturbationElement*             freeEnergyPerturbationElementPtr,
                    bool                                       hasReadEkinState);

    //! Build the force element - can be normal forces or shell / flex constraints
    std::unique_ptr<ISimulatorElement>
    buildForces(SignallerBuilder<NeighborSearchSignaller>* neighborSearchSignallerBuilder,
                SignallerBuilder<EnergySignaller>*         energySignallerBuilder,
                StatePropagatorData*                       statePropagatorDataPtr,
                EnergyElement*                             energyElementPtr,
                FreeEnergyPerturbationElement*             freeEnergyPerturbationElement);

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
     *       repeatedly. While preliminary measures do not seem to indicate
     *       performance problems, a pool allocator would be an obvious
     *       optimization possibility, as would reusing the task queue if
     *       the task queue is periodic around the rebuild point (i.e. the
     *       task queue is identical between rebuilds).
     */
    void populateTaskQueue();

    //! Check for disabled functionality (during construction time)
    void checkInputForDisabledFunctionality();

    //! The run queue
    std::queue<SimulatorRunFunctionPtr> taskQueue_;

    /* Note that the Simulator is owning the signallers and elements.
     * The ownership list and the call list are kept separate, however,
     * to allow to have elements more than once in the call lists -
     * either as signaller AND element (such as the TrajectoryElement),
     * or to have an element twice in the scheduling sequence (currently
     * not used).
     *
     * For the elements, the setup and teardown is applied on the
     * elementsOwnershipList_, to ensure that it is called only once per
     * element. For the signallers, the setup is applied on the
     * signallerCallList_ - this makes sure that both the elementSetup()
     * and signallerSetup() of an object being both an element and a
     * signaller is called. It is also not expected to run have a signaller
     * more than once in the signallerCallList_, so we don't have to worry
     * about calling the setup method twice. Consequently, this means that
     * objects being both a signaller and an element should be stored in
     * the elementsOwnershipList_.
     */
    //! List of signalers (ownership)
    std::vector<std::unique_ptr<ISignaller>> signallersOwnershipList_;
    //! List of signalers (calling sequence)
    std::vector<compat::not_null<ISignaller*>> signallerCallList_;
    //! List of schedulerElements (ownership)
    std::vector<std::unique_ptr<ISimulatorElement>> elementsOwnershipList_;
    //! List of schedulerElements (calling sequence)
    std::vector<compat::not_null<ISimulatorElement*>> elementCallList_;

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
    //! Compute globals communication period
    int nstglobalcomm_;

    //! The topology
    std::unique_ptr<TopologyHolder> topologyHolder_;

    //! The current step
    Step step_ = -1;

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
};

//! Constructor implementation (here to avoid template-related linker problems)
template<typename... Args>
ModularSimulator::ModularSimulator(Args&&... args) : ISimulator(std::forward<Args>(args)...)
{
    nstglobalcomm_ = computeGlobalCommunicationPeriod(mdlog, inputrec, cr);
    signalHelper_  = std::make_unique<SignalHelper>();
    checkInputForDisabledFunctionality();
}

//! \cond
template<typename T, typename U>
void ModularSimulator::addToCallList(U* element, std::vector<compat::not_null<T*>>& callList)
{
    if (element)
    {
        callList.emplace_back(element);
    }
}

template<typename T, typename U>
void ModularSimulator::addToCallList(gmx::compat::not_null<U*>          element,
                                     std::vector<compat::not_null<T*>>& callList)
{
    callList.emplace_back(element);
}

template<typename T, typename U>
void ModularSimulator::addToCallList(std::unique_ptr<U>& element, std::vector<compat::not_null<T*>>& callList)
{
    if (element)
    {
        callList.emplace_back(compat::make_not_null(element.get()));
    }
}

template<typename T, typename U>
void ModularSimulator::addToCallListAndMove(std::unique_ptr<U>                 element,
                                            std::vector<compat::not_null<T*>>& callList,
                                            std::vector<std::unique_ptr<T>>&   elementList)
{
    if (element)
    {
        callList.emplace_back(compat::make_not_null(element.get()));
        elementList.emplace_back(std::move(element));
    }
}
//! \endcond

/*!
 * \brief Whether or not to use the ModularSimulator
 *
 * GMX_DISABLE_MODULAR_SIMULATOR environment variable allows to disable modular simulator for
 * all uses.
 *
 * See ModularSimulator::isInputCompatible() for function signature.
 *
 * \ingroup module_modularsimulator
 */
template<typename... Ts>
auto checkUseModularSimulator(Ts&&... args)
        -> decltype(ModularSimulator::isInputCompatible(std::forward<Ts>(args)...))
{
    return ModularSimulator::isInputCompatible(std::forward<Ts>(args)...)
           && getenv("GMX_DISABLE_MODULAR_SIMULATOR") == nullptr;
}

} // namespace gmx

#endif // GROMACS_MODULARSIMULATOR_MODULARSIMULATOR_H
