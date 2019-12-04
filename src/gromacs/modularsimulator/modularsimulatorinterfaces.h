/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \defgroup module_modularsimulator The modular simulator
 * \ingroup group_mdrun
 * \brief
 * The modular simulator improves extensibility, adds Monte Carlo capabilities,
 * promotes data locality and communication via interfaces, supports
 * multi-stepping integrators, and paves the way for some task parallelism.
 *
 * For more information, see page_modularsimulator
 * \todo Can we link to `docs/doxygen/lib/modularsimulator.md`?
 *
 * \author Pascal Merz <pascal.merz@me.com>
 */
/*! \libinternal \file
 * \brief
 * Declares the main interfaces used by the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */
#ifndef GMX_MODULARSIMULATOR_MODULARSIMULATORINTERFACES_H
#define GMX_MODULARSIMULATOR_MODULARSIMULATORINTERFACES_H

#include <functional>
#include <memory>

#include "gromacs/utility/basedefinitions.h"

struct gmx_localtop_t;
struct gmx_mdoutf;
class t_state;

namespace gmx
{
template<class Signaller>
class SignallerBuilder;
class NeighborSearchSignaller;
class LastStepSignaller;
class LoggingSignaller;
class EnergySignaller;

//! \addtogroup module_modularsimulator
//! \{

// This will make signatures more legible
//! Step number
using Step = int64_t;
//! Simulation time
using Time = double;

//! The function type that can be scheduled to be run during the simulator run
typedef std::function<void()> SimulatorRunFunction;
//! Pointer to the function type that can be scheduled to be run during the simulator run
typedef std::unique_ptr<SimulatorRunFunction> SimulatorRunFunctionPtr;

//! The function type that allows to register run functions
typedef std::function<void(SimulatorRunFunctionPtr)> RegisterRunFunction;
//! Pointer to the function type that allows to register run functions
typedef std::unique_ptr<RegisterRunFunction> RegisterRunFunctionPtr;

/*! \libinternal
 * \brief The general interface for elements of the modular simulator
 *
 * Setup and teardown are run once at the beginning of the simulation
 * (after all elements were set up, before the first step) and at the
 * end of the simulation (after the last step, before elements are
 * destructed), respectively. `registerRun` is periodically called
 * during the run, at which point elements can decide to register one
 * or more run functions to be run at a specific time / step. Registering
 * more than one function is especially valuable for collective elements
 * consisting of more than one single element.
 */
class ISimulatorElement
{
public:
    /*! \brief Query whether element wants to run at step / time
     *
     * Element can register one or more functions to be run at that step through
     * the registration pointer.
     */
    virtual void scheduleTask(Step, Time, const RegisterRunFunctionPtr&) = 0;
    //! Method guaranteed to be called after construction, before simulator run
    virtual void elementSetup() = 0;
    //! Method guaranteed to be called after simulator run, before deconstruction
    virtual void elementTeardown() = 0;
    //! Standard virtual destructor
    virtual ~ISimulatorElement() = default;
};

/*! \libinternal
 * \brief The general Signaller interface
 *
 * Signallers are run at the beginning of Simulator steps, informing
 * their clients about the upcoming step. This allows clients to
 * decide if they need to be activated at this time step, and what
 * functions they will have to run. Examples for signallers
 * include the neighbor-search signaller (informs its clients when a
 * neighbor-searching step is about to happen) or the logging
 * signaller (informs its clients when a logging step is about to
 * happen).
 *
 * We expect all signallers to accept registration from clients, but
 * different signallers might handle that differently, so we don't
 * include this in the interface.
 */
class ISignaller
{
public:
    //! Function run before every step of scheduling
    virtual void signal(Step, Time) = 0;
    //! Method guaranteed to be called after construction, before simulator run
    virtual void signallerSetup() = 0;
    //! Standard virtual destructor
    virtual ~ISignaller() = default;
};

//! The function type that can be registered to signallers for callback
typedef std::function<void(Step, Time)> SignallerCallback;
//! Pointer to the function type that can be registered to signallers for callback
typedef std::unique_ptr<SignallerCallback> SignallerCallbackPtr;

/*! \libinternal
 * \brief Interface for clients of the NeighborSearchSignaller
 *
 * Defining registerNSCallback allows clients to register an arbitrary callback
 * for notification by the signaller.
 */
class INeighborSearchSignallerClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow builder of NeighborSearchSignaller to ask for callback registration
    friend class SignallerBuilder<NeighborSearchSignaller>;
    //! @endcond
    //! Standard virtual destructor
    virtual ~INeighborSearchSignallerClient() = default;

protected:
    //! Return callback to NeighborSearchSignaller
    virtual SignallerCallbackPtr registerNSCallback() = 0;
};

/*! \libinternal
 * \brief Interface for clients of the LastStepSignaller
 *
 * Defining registerLastStepCallback allows clients to register an arbitrary callback
 * for notification by the signaller.
 */
class ILastStepSignallerClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow builder of LastStepSignaller to ask for callback registration
    friend class SignallerBuilder<LastStepSignaller>;
    //! @endcond
    //! Standard virtual destructor
    virtual ~ILastStepSignallerClient() = default;

protected:
    //! Return callback to LastStepSignaller
    virtual SignallerCallbackPtr registerLastStepCallback() = 0;
};

/*! \libinternal
 * \brief Interface for clients of the LoggingSignaller
 *
 * Defining registerLoggingCallback allows clients to register an arbitrary callback
 * for notification by the signaller.
 */
class ILoggingSignallerClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow builder of LoggingSignaller to ask for callback registration
    friend class SignallerBuilder<LoggingSignaller>;
    //! @endcond
    //! Standard virtual destructor
    virtual ~ILoggingSignallerClient() = default;

protected:
    //! Return callback to LoggingSignaller
    virtual SignallerCallbackPtr registerLoggingCallback() = 0;
};

//! The energy events signalled by the EnergySignaller
enum class EnergySignallerEvent
{
    EnergyCalculationStep,
    VirialCalculationStep,
    FreeEnergyCalculationStep
};

/*! \libinternal
 * \brief Interface for clients of the EnergySignaller
 *
 * Defining registerEnergyCallback allows clients to register an arbitrary callback
 * for notification by the signaller for every EnergySignallerEvent separately.
 */
class IEnergySignallerClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow builder of EnergySignaller to ask for callback registration
    friend class SignallerBuilder<EnergySignaller>;
    //! @endcond
    //! Standard virtual destructor
    virtual ~IEnergySignallerClient() = default;

protected:
    //! Return callback to EnergySignaller
    virtual SignallerCallbackPtr registerEnergyCallback(EnergySignallerEvent) = 0;
};

//! The trajectory writing events
enum class TrajectoryEvent
{
    StateWritingStep,
    EnergyWritingStep
};

/*! \libinternal
 * \brief Interface for signaller clients of the TrajectoryElement
 *
 * Defining registerTrajectorySignallerCallback allows clients to register an arbitrary
 * callback for notification by the signaller for every TrajectoryEvent separately.
 */
class ITrajectorySignallerClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow builder of TrajectoryElement to ask for callback registration
    friend class TrajectoryElementBuilder;
    //! @endcond
    //! Standard virtual destructor
    virtual ~ITrajectorySignallerClient() = default;

protected:
    //! Return callback to TrajectoryElement
    virtual SignallerCallbackPtr registerTrajectorySignallerCallback(TrajectoryEvent) = 0;
};

/* Trajectory writing clients are handed a pointer to the output file handler,
 * allowing them to write their own trajectory contribution.
 *
 * As trajectory writing and log writing cannot currently be separated for the
 * energy, clients also get informed whether this is a trajectory-writing step
 * and / or a log-writing step.
 */
//! Function type for trajectory writing clients
typedef std::function<void(gmx_mdoutf*, Step, Time, bool, bool)> ITrajectoryWriterCallback;
//! Pointer to the function type for trajectory writing clients
typedef std::unique_ptr<ITrajectoryWriterCallback> ITrajectoryWriterCallbackPtr;

/*! \libinternal
 * \brief Interface for writer clients of the TrajectoryElement
 *
 * Defining registerTrajectoryWriterCallback allows clients to register an arbitrary
 * callback called by the TrajectoryElement when trajectory writing happens.
 *
 * Setup and teardown methods allow clients to perform tasks with a valid output pointer.
 */
class ITrajectoryWriterClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow TrajectoryElement to ask for callback registration
    friend class TrajectoryElement;
    //! @endcond
    //! Standard virtual destructor
    virtual ~ITrajectoryWriterClient() = default;

protected:
    //! Setup method with valid output pointer.
    virtual void trajectoryWriterSetup(gmx_mdoutf* outf) = 0;
    //! Teardown method with valid output pointer.
    virtual void trajectoryWriterTeardown(gmx_mdoutf* outf) = 0;

    //! Return callback to TrajectoryElement
    virtual ITrajectoryWriterCallbackPtr registerTrajectoryWriterCallback(TrajectoryEvent) = 0;
};

/*! \libinternal
 * \brief Client requiring read access to the local topology
 *
 */
class ITopologyHolderClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow TopologyHolder to set new topology
    friend class TopologyHolder;
    //! @endcond
    //! Standard virtual destructor
    virtual ~ITopologyHolderClient() = default;

protected:
    //! Pass pointer to new local topology
    virtual void setTopology(const gmx_localtop_t*) = 0;
};

/*! \libinternal
 * \brief Client that needs to store data during checkpointing
 *
 * The current checkpointing helper uses the legacy t_state object to collect
 * the data to be checkpointed. Clients get queried for their contributions
 * using pointers to t_state objects.
 * \todo Add checkpoint reading
 * \todo Evolve this to a model in which the checkpoint helper passes a file
 *       pointer rather than a t_state object, and the clients are responsible
 *       to read / write.
 */
class ICheckpointHelperClient
{
public:
    //! @cond
    // (doxygen doesn't like these...)
    //! Allow CheckpointHelper to interact
    friend class CheckpointHelper;
    //! @endcond
    //! Standard virtual destructor
    virtual ~ICheckpointHelperClient() = default;

protected:
    //! Write checkpoint
    virtual void writeCheckpoint(t_state* localState, t_state* globalState) = 0;
};
//! /}
} // namespace gmx

#endif // GMX_MODULARSIMULATOR_MODULARSIMULATORINTERFACES_H
