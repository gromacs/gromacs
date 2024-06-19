/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Declares the signallers for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_SIGNALLERS_H
#define GMX_MODULARSIMULATOR_SIGNALLERS_H

#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "gromacs/compat/pointers.h"

#include "modularsimulatorinterfaces.h"

namespace gmx
{
class StopHandler;
class TrajectoryElement;
enum class StartingBehavior;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Builder for signallers
 *
 * This builder allows clients to register, and then builds the signaller
 * passing on the list of clients.
 *
 * \tparam Signaller  The signaller to be built
 */
template<typename Signaller>
class SignallerBuilder final
{
public:
    //! Allows clients to register to the signaller
    void registerSignallerClient(typename Signaller::Client* client);

    //! Build the signaller
    template<typename... Args>
    std::unique_ptr<Signaller> build(Args&&... args);

private:
    //! List of signaller clients
    std::vector<typename Signaller::Client*> signallerClients_;
    //! The state of the builder
    ModularSimulatorBuilderState state_ = ModularSimulatorBuilderState::AcceptingClientRegistrations;

    //! Helper function to get the callbacks from the clients
    template<typename... Args>
    std::vector<SignallerCallback> buildCallbackVector(Args&&... args);

    /*! \brief Get a callback from a single client
     *
     * This is in a separate function, as the exact call depends on the
     * specific signaller / client.
     */
    template<typename... Args>
    std::optional<SignallerCallback> getSignallerCallback(typename Signaller::Client* client,
                                                          Args&&... args);
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element signalling a neighbor search step
 *
 * This element informs its clients via callbacks
 * when a neighbor-searching step is happening.
 */
class NeighborSearchSignaller final : public ISignaller
{
public:
    /*! \brief Run the signaller at a specific step / time
     *
     * Informs callbacks if step % nstlist_ == 0
     *
     * \param step  The current time step
     * \param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Do nothing at setup time
    void setup() override{};

    //! Allow builder to construct
    friend class SignallerBuilder<NeighborSearchSignaller>;
    //! Define client type
    typedef INeighborSearchSignallerClient Client;

private:
    /*! \brief Constructor
     *
     * \param callbacks  A vector of pointers to callbacks
     * \param nstlist    The frequency at which neighbor search is performed
     * \param initStep   The first step of the simulation
     * \param initTime   The start time of the simulation
     */
    NeighborSearchSignaller(std::vector<SignallerCallback> callbacks, Step nstlist, Step initStep, Time initTime);

    //! Client callbacks
    std::vector<SignallerCallback> callbacks_;

    //! The NS frequency
    const Step nstlist_;
    //! The initial step of the simulation
    const Step initStep_;
    //! The initial time of the simulation
    const Time initTime_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element signalling the last step
 *
 * This element informs its clients via callbacks
 * when the last step is happening.
 */
class LastStepSignaller final : public ISignaller, public INeighborSearchSignallerClient
{
public:
    /*! \brief Run the signaller at a specific step / time
     *
     * Informs callbacks if this is the last step
     *
     * \param step  The current time step
     * \param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Check that necessary registration was done
    void setup() override;

    //! Allow builder to construct
    friend class SignallerBuilder<LastStepSignaller>;
    //! Define client type
    typedef ILastStepSignallerClient Client;

private:
    /*! \brief Constructor
     *
     * \param callbacks    A vector of pointers to callbacks
     * \param nsteps       The total number of steps for the simulation
     * \param initStep     The first step of the simulation
     * \param stopHandler  A pointer to the stop handler (LastStepSignaller takes ownership)
     */
    LastStepSignaller(std::vector<SignallerCallback> callbacks, Step nsteps, Step initStep, StopHandler* stopHandler);

    //! Client callbacks
    std::vector<SignallerCallback> callbacks_;

    //! The last step of the simulation
    const Step stopStep_;
    //! Whether we signalled last step due to stop condition
    bool signalledStopCondition_;
    //! A pointer to the stop handler communicating signal and time-related stops
    StopHandler* stopHandler_;

    //! INeighborSearchSignallerClient implementation
    std::optional<SignallerCallback> registerNSCallback() override;
    //! The next NS step (notified by NS signaller)
    Step nextNSStep_;
    //! Whether we registered to the NS signaller
    bool nsStepRegistrationDone_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element signalling a logging step
 *
 * This element informs its clients via callbacks
 * when a logging step is happening.
 */
class LoggingSignaller final : public ISignaller, public ILastStepSignallerClient
{
public:
    /*! \brief Run the signaller at a specific step / time
     *
     * Informs callbacks if step % nstlog_ == 0
     *
     * \param step  The current time step
     * \param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Check that necessary registration was done
    void setup() override;

    //! Allow builder to construct
    friend class SignallerBuilder<LoggingSignaller>;
    //! Define client type
    typedef ILoggingSignallerClient Client;

private:
    /*! \brief Constructor
     *
     * \param callbacks         A vector of pointers to callbacks
     * \param nstlog            The logging frequency
     * \param initStep          The first step of the simulation
     * \param startingBehavior  Whether this is a new simulation or restarting from checkpoint
     */
    LoggingSignaller(std::vector<SignallerCallback> callbacks,
                     Step                           nstlog,
                     Step                           initStep,
                     StartingBehavior               startingBehavior);

    //! Client callbacks
    std::vector<SignallerCallback> callbacks_;

    //! The logging frequency
    const Step nstlog_;
    //! The initial step of the simulation
    const Step initStep_;
    //! How we are starting the simulation
    const StartingBehavior startingBehavior_;

    //! ILastStepSignallerClient implementation
    std::optional<SignallerCallback> registerLastStepCallback() override;
    //! The last step (notified by signaller)
    Step lastStep_;
    //! Whether we registered to the last step signaller
    bool lastStepRegistrationDone_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element signalling trajectory writing
 *
 * During signalling phase, it checks whether the current step is a writing
 * step for either the energy or the state (position, velocity, forces)
 * trajectory. It then notifies the signaller clients of the upcoming step.
 *
 * The TrajectorySignaller works in close collaboration with the TrajectoryElement
 * which does the actual trajectory writing during the simulation step.
 */
class TrajectorySignaller final : public ISignaller, public ILastStepSignallerClient
{
public:
    /*! \brief Prepare signaller
     *
     * Check that necessary registration was done
     */
    void setup() override;

    /*! \brief Run the signaller at a specific step / time
     *
     * Informs clients when energy or state will be written.
     *
     * \param step           The current time step
     * \param time           The current time
     */
    void signal(Step step, Time time) override;

    //! Allow builder to construct
    friend class SignallerBuilder<TrajectorySignaller>;
    //! Define client type
    typedef ITrajectorySignallerClient Client;

private:
    //! Constructor
    TrajectorySignaller(std::vector<SignallerCallback> signalEnergyCallbacks,
                        std::vector<SignallerCallback> signalStateCallbacks,
                        int                            nstxout,
                        int                            nstvout,
                        int                            nstfout,
                        int                            nstxoutCompressed,
                        int                            tngBoxOut,
                        int                            tngLambdaOut,
                        int                            tngBoxOutCompressed,
                        int                            tngLambdaOutCompressed,
                        int                            nstenergy);

    //! Output frequencies
    //! {
    const int nstxout_;
    const int nstvout_;
    const int nstfout_;
    const int nstxoutCompressed_;
    const int tngBoxOut_;
    const int tngLambdaOut_;
    const int tngBoxOutCompressed_;
    const int tngLambdaOutCompressed_;
    const int nstenergy_;
    //! }

    //! Callbacks to signal events
    //! {
    std::vector<SignallerCallback> signalEnergyCallbacks_;
    std::vector<SignallerCallback> signalStateCallbacks_;
    //! }

    /*
     * Last step client
     */
    Step lastStep_;
    bool lastStepRegistrationDone_;
    //! ILastStepSignallerClient implementation
    std::optional<SignallerCallback> registerLastStepCallback() override;
};

//! When we calculate virial
enum class EnergySignallerVirialMode
{
    Off,           //!< No specific virial calculation - calculate when energy is calculated
    OnStep,        //!< Calculate on virial frequency steps
    OnStepAndNext, //!< Calculate on virial frequency steps and on step after
    Count          //!< The number of entries
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element signalling energy related special steps
 *
 * This element informs its clients via callbacks
 * of the following events:
 *   - energy calculation step
 *   - virial calculation step
 *   - free energy calculation step
 */
class EnergySignaller final : public ISignaller, public ITrajectorySignallerClient, public ILoggingSignallerClient
{
public:
    /*! \brief Run the signaller at a specific step / time
     *
     * Informs callbacks of energy / virial / free energy special steps
     *
     * \param step  The current time step
     * \param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Check that necessary registration was done
    void setup() override;

    //! Allow builder to construct
    friend class SignallerBuilder<EnergySignaller>;
    //! Define client type
    typedef IEnergySignallerClient Client;

private:
    /*! \brief Constructor
     *
     * \param calculateEnergyCallbacks      A vector of pointers to callbacks (energy steps)
     * \param calculateVirialCallbacks      A vector of pointers to callbacks (virial steps)
     * \param calculateFreeEnergyCallbacks  A vector of pointers to callbacks (free energy steps)
     * \param nstcalcenergy                 The energy calculation frequency
     * \param nstcalcfreeenergy             The free energy calculation frequency
     * \param nstcalcvirial                 The free energy calculation frequency
     * \param virialMode                    Indicates which steps will calculate virial
     */
    EnergySignaller(std::vector<SignallerCallback> calculateEnergyCallbacks,
                    std::vector<SignallerCallback> calculateVirialCallbacks,
                    std::vector<SignallerCallback> calculateFreeEnergyCallbacks,
                    int                            nstcalcenergy,
                    int                            nstcalcfreeenergy,
                    int                            nstcalcvirial,
                    EnergySignallerVirialMode      virialMode);

    //! Client callbacks
    //! {
    std::vector<SignallerCallback> calculateEnergyCallbacks_;
    std::vector<SignallerCallback> calculateVirialCallbacks_;
    std::vector<SignallerCallback> calculateFreeEnergyCallbacks_;
    //! }

    //! The energy calculation frequency
    const int nstcalcenergy_;
    //! The free energy calculation frequency
    const int nstcalcfreeenergy_;
    //! The virial calculation frequency
    const int nstcalcvirial_;
    //! The virial calculation mode
    const EnergySignallerVirialMode virialMode_;

    //! ITrajectorySignallerClient implementation
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! The energy writing step (notified by signaller)
    Step energyWritingStep_;
    //! Whether we registered to the trajectory signaller
    bool trajectoryRegistrationDone_;

    //! ILoggingSignallerClient implementation
    std::optional<SignallerCallback> registerLoggingCallback() override;
    //! The next logging step (notified by signaller)
    Step loggingStep_;
    //! Whether we registered to the logging signaller
    bool loggingRegistrationDone_;
};

//! Allows clients to register to the signaller
template<class Signaller>
void SignallerBuilder<Signaller>::registerSignallerClient(typename Signaller::Client* client)
{
    if (client)
    {
        if (state_ == ModularSimulatorBuilderState::NotAcceptingClientRegistrations)
        {
            throw SimulationAlgorithmSetupError(
                    "Tried to register to signaller after it was built.");
        }
        signallerClients_.emplace_back(client);
    }
}

/*! \brief Build the signaller
 *
 * General version - for NeighborSearchSignaller, LastStepSignaller, LoggingSignaller
 */
template<class Signaller>
template<typename... Args>
std::unique_ptr<Signaller> SignallerBuilder<Signaller>::build(Args&&... args)
{
    state_         = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
    auto callbacks = buildCallbackVector();
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<Signaller>(new Signaller(std::move(callbacks), std::forward<Args>(args)...));
}

/*! \brief Build the signaller
 *
 * Specialized version - TrajectorySignaller has a different build process
 */
template<>
template<typename... Args>
std::unique_ptr<TrajectorySignaller> SignallerBuilder<TrajectorySignaller>::build(Args&&... args)
{
    state_                     = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
    auto signalEnergyCallbacks = buildCallbackVector(TrajectoryEvent::EnergyWritingStep);
    auto signalStateCallbacks  = buildCallbackVector(TrajectoryEvent::StateWritingStep);
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<TrajectorySignaller>(new TrajectorySignaller(
            std::move(signalEnergyCallbacks), std::move(signalStateCallbacks), std::forward<Args>(args)...));
}

/*! \brief Build the signaller
 *
 * Specialized version - EnergySignaller has a significantly different build process
 */
template<>
template<typename... Args>
std::unique_ptr<EnergySignaller> SignallerBuilder<EnergySignaller>::build(Args&&... args)
{
    state_                        = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
    auto calculateEnergyCallbacks = buildCallbackVector(EnergySignallerEvent::EnergyCalculationStep);
    auto calculateVirialCallbacks = buildCallbackVector(EnergySignallerEvent::VirialCalculationStep);
    auto calculateFreeEnergyCallbacks =
            buildCallbackVector(EnergySignallerEvent::FreeEnergyCalculationStep);
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<EnergySignaller>(new EnergySignaller(std::move(calculateEnergyCallbacks),
                                                                std::move(calculateVirialCallbacks),
                                                                std::move(calculateFreeEnergyCallbacks),
                                                                std::forward<Args>(args)...));
}

//! Helper function to get the callbacks from the clients
template<typename Signaller>
template<typename... Args>
std::vector<SignallerCallback> SignallerBuilder<Signaller>::buildCallbackVector(Args&&... args)
{
    std::vector<SignallerCallback> callbacks;
    // Allow clients to register their callbacks
    for (auto& client : signallerClients_)
    {
        if (auto callback = getSignallerCallback(client, std::forward<Args>(args)...)) // don't register nullptr
        {
            callbacks.emplace_back(std::move(*callback));
        }
    }
    return callbacks;
}

//! Get a callback from a single client - NeighborSearchSignaller
template<>
template<typename... Args>
std::optional<SignallerCallback> SignallerBuilder<NeighborSearchSignaller>::getSignallerCallback(
        typename NeighborSearchSignaller::Client* client,
        Args&&... args)
{
    return client->registerNSCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - LastStepSignaller
template<>
template<typename... Args>
std::optional<SignallerCallback>
SignallerBuilder<LastStepSignaller>::getSignallerCallback(typename LastStepSignaller::Client* client,
                                                          Args&&... args)
{
    return client->registerLastStepCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - LoggingSignaller
template<>
template<typename... Args>
std::optional<SignallerCallback>
SignallerBuilder<LoggingSignaller>::getSignallerCallback(typename LoggingSignaller::Client* client,
                                                         Args&&... args)
{
    return client->registerLoggingCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - TrajectorySignaller
template<>
template<typename... Args>
std::optional<SignallerCallback>
SignallerBuilder<TrajectorySignaller>::getSignallerCallback(typename TrajectorySignaller::Client* client,
                                                            Args&&... args)
{
    return client->registerTrajectorySignallerCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - EnergySignaller
template<>
template<typename... Args>
std::optional<SignallerCallback>
SignallerBuilder<EnergySignaller>::getSignallerCallback(typename EnergySignaller::Client* client,
                                                        Args&&... args)
{
    return client->registerEnergyCallback(std::forward<Args>(args)...);
}

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_SIGNALLERS_H
