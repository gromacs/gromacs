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
/*! \libinternal \file
 * \brief Declares the signallers for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_SIGNALLERS_H
#define GMX_MODULARSIMULATOR_SIGNALLERS_H

#include <vector>

#include "gromacs/compat/pointers.h"

#include "modularsimulatorinterfaces.h"

namespace gmx
{
class StopHandler;
class TrajectoryElement;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Builder for signallers
 *
 * This builder allows clients to register, and then builds the signaller
 * passing on the list of clients.
 *
 * @tparam Signaller  The signaller to be built
 */
template<typename Signaller>
class SignallerBuilder final
{
public:
    //! Allows clients to register to the signaller
    void registerSignallerClient(compat::not_null<typename Signaller::Client*> client);

    //! Build the signaller
    template<typename... Args>
    std::unique_ptr<Signaller> build(Args&&... args);

private:
    //! List of signaller clients
    std::vector<typename Signaller::Client*> signallerClients_;

    //! Helper function to get the callbacks from the clients
    template<typename... Args>
    std::vector<SignallerCallbackPtr> buildCallbackVector(Args&&... args);

    /*! \brief Get a callback from a single client
     *
     * This is in a separate function, as the exact call depends on the
     * specific signaller / client.
     */
    template<typename... Args>
    SignallerCallbackPtr getSignallerCallback(typename Signaller::Client* client, Args&&... args);
};

/*! \libinternal
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
     * @param step  The current time step
     * @param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Do nothing at setup time
    void signallerSetup() override{};

    //! Allow builder to construct
    friend class SignallerBuilder<NeighborSearchSignaller>;
    //! Define client type
    typedef INeighborSearchSignallerClient Client;

private:
    /*! \brief Constructor
     *
     * @param callbacks  A vector of pointers to callbacks
     * @param nstlist    The frequency at which neighbor search is performed
     * @param initStep   The first step of the simulation
     * @param initTime   The start time of the simulation
     */
    NeighborSearchSignaller(std::vector<SignallerCallbackPtr> callbacks, Step nstlist, Step initStep, Time initTime);

    //! Client callbacks
    std::vector<SignallerCallbackPtr> callbacks_;

    //! The NS frequency
    const Step nstlist_;
    //! The initial step of the simulation
    const Step initStep_;
    //! The initial time of the simulation
    const Time initTime_;
};

/*! \libinternal
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
     * @param step  The current time step
     * @param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Check that necessary registration was done
    void signallerSetup() override;

    //! Allow builder to construct
    friend class SignallerBuilder<LastStepSignaller>;
    //! Define client type
    typedef ILastStepSignallerClient Client;

private:
    /*! \brief Constructor
     *
     * @param callbacks    A vector of pointers to callbacks
     * @param nsteps       The total number of steps for the simulation
     * @param initStep     The first step of the simulation
     * @param stopHandler  A pointer to the stop handler (LastStepSignaller takes ownership)
     */
    LastStepSignaller(std::vector<SignallerCallbackPtr> callbacks,
                      Step                              nsteps,
                      Step                              initStep,
                      StopHandler*                      stopHandler);

    //! Client callbacks
    std::vector<SignallerCallbackPtr> callbacks_;

    //! The last step of the simulation
    const Step stopStep_;
    //! Whether we signalled last step due to stop condition
    bool signalledStopCondition_;
    //! A pointer to the stop handler communicating signal and time-related stops
    StopHandler* stopHandler_;

    //! INeighborSearchSignallerClient implementation
    SignallerCallbackPtr registerNSCallback() override;
    //! The next NS step (notified by NS signaller)
    Step nextNSStep_;
    //! Whether we registered to the NS signaller
    bool nsStepRegistrationDone_;
};

/*! \libinternal
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
     * @param step  The current time step
     * @param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Check that necessary registration was done
    void signallerSetup() override;

    //! Allow builder to construct
    friend class SignallerBuilder<LoggingSignaller>;
    //! Define client type
    typedef ILoggingSignallerClient Client;

private:
    /*! \brief Constructor
     *
     * @param callbacks  A vector of pointers to callbacks
     * @param nstlog     The logging frequency
     * @param initStep   The first step of the simulation
     * @param initTime   The start time of the simulation
     */
    LoggingSignaller(std::vector<SignallerCallbackPtr> callbacks, Step nstlog, Step initStep, Time initTime);

    //! Client callbacks
    std::vector<SignallerCallbackPtr> callbacks_;

    //! The logging frequency
    const Step nstlog_;
    //! The initial step of the simulation
    const Step initStep_;
    //! The initial time of the simulation
    const Time initTime_;

    //! ILastStepSignallerClient implementation
    SignallerCallbackPtr registerLastStepCallback() override;
    //! The last step (notified by signaller)
    Step lastStep_;
    //! Whether we registered to the last step signaller
    bool lastStepRegistrationDone_;
};

/*! \libinternal
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
     * @param step  The current time step
     * @param time  The current time
     */
    void signal(Step step, Time time) override;

    //! Check that necessary registration was done
    void signallerSetup() override;

    //! Allow builder to construct
    friend class SignallerBuilder<EnergySignaller>;
    //! Define client type
    typedef IEnergySignallerClient Client;

private:
    /*! \brief Constructor
     *
     * @param calculateEnergyCallbacks      A vector of pointers to callbacks (energy steps)
     * @param calculateVirialCallbacks      A vector of pointers to callbacks (virial steps)
     * @param calculateFreeEnergyCallbacks  A vector of pointers to callbacks (free energy steps)
     * @param nstcalcenergy                 The energy calculation frequency
     * @param nstcalcfreeenergy             The free energy calculation frequency
     * @param nstcalcvirial                 The free energy calculation frequency
     */
    EnergySignaller(std::vector<SignallerCallbackPtr> calculateEnergyCallbacks,
                    std::vector<SignallerCallbackPtr> calculateVirialCallbacks,
                    std::vector<SignallerCallbackPtr> calculateFreeEnergyCallbacks,
                    int                               nstcalcenergy,
                    int                               nstcalcfreeenergy,
                    int                               nstcalcvirial);

    //! Client callbacks
    //! {
    std::vector<SignallerCallbackPtr> calculateEnergyCallbacks_;
    std::vector<SignallerCallbackPtr> calculateVirialCallbacks_;
    std::vector<SignallerCallbackPtr> calculateFreeEnergyCallbacks_;
    //! }

    //! The energy calculation frequency
    const int nstcalcenergy_;
    //! The free energy calculation frequency
    const int nstcalcfreeenergy_;
    //! The virial calculation frequency
    const int nstcalcvirial_;

    //! ITrajectorySignallerClient implementation
    SignallerCallbackPtr registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! The energy writing step (notified by signaller)
    Step energyWritingStep_;
    //! Whether we registered to the trajectory signaller
    bool trajectoryRegistrationDone_;

    //! ILoggingSignallerClient implementation
    SignallerCallbackPtr registerLoggingCallback() override;
    //! The next logging step (notified by signaller)
    Step loggingStep_;
    //! Whether we registered to the logging signaller
    bool loggingRegistrationDone_;
};

//! Allows clients to register to the signaller
template<class Signaller>
void SignallerBuilder<Signaller>::registerSignallerClient(compat::not_null<typename Signaller::Client*> client)
{
    signallerClients_.emplace_back(client);
}

/*! \brief Build the signaller
 *
 * General version - for NeighborSearchSignaller, LastStepSignaller, LoggingSignaller
 */
template<class Signaller>
template<typename... Args>
std::unique_ptr<Signaller> SignallerBuilder<Signaller>::build(Args&&... args)
{
    auto callbacks = buildCallbackVector();
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<Signaller>(new Signaller(std::move(callbacks), std::forward<Args>(args)...));
}

/*! \brief Build the signaller
 *
 * Specialized version - EnergySignaller has a significantly different build process
 */
template<>
template<typename... Args>
std::unique_ptr<EnergySignaller> SignallerBuilder<EnergySignaller>::build(Args&&... args)
{
    auto calculateEnergyCallbacks = buildCallbackVector(EnergySignallerEvent::EnergyCalculationStep);
    auto calculateVirialCallbacks = buildCallbackVector(EnergySignallerEvent::VirialCalculationStep);
    auto calculateFreeEnergyCallbacks =
            buildCallbackVector(EnergySignallerEvent::FreeEnergyCalculationStep);
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<EnergySignaller>(new EnergySignaller(
            std::move(calculateEnergyCallbacks), std::move(calculateVirialCallbacks),
            std::move(calculateFreeEnergyCallbacks), std::forward<Args>(args)...));
}

//! Helper function to get the callbacks from the clients
template<typename Signaller>
template<typename... Args>
std::vector<SignallerCallbackPtr> SignallerBuilder<Signaller>::buildCallbackVector(Args&&... args)
{
    std::vector<SignallerCallbackPtr> callbacks;
    // Allow clients to register their callbacks
    for (auto& client : signallerClients_)
    {
        if (auto callback = getSignallerCallback(client, std::forward<Args>(args)...)) // don't register nullptr
        {
            callbacks.emplace_back(std::move(callback));
        }
    }
    return callbacks;
}

//! Get a callback from a single client - NeighborSearchSignaller
template<>
template<typename... Args>
SignallerCallbackPtr SignallerBuilder<NeighborSearchSignaller>::getSignallerCallback(
        typename NeighborSearchSignaller::Client* client,
        Args&&... args)
{
    return client->registerNSCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - LastStepSignaller
template<>
template<typename... Args>
SignallerCallbackPtr
SignallerBuilder<LastStepSignaller>::getSignallerCallback(typename LastStepSignaller::Client* client,
                                                          Args&&... args)
{
    return client->registerLastStepCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - LoggingSignaller
template<>
template<typename... Args>
SignallerCallbackPtr
SignallerBuilder<LoggingSignaller>::getSignallerCallback(typename LoggingSignaller::Client* client,
                                                         Args&&... args)
{
    return client->registerLoggingCallback(std::forward<Args>(args)...);
}

//! Get a callback from a single client - EnergySignaller
template<>
template<typename... Args>
SignallerCallbackPtr SignallerBuilder<EnergySignaller>::getSignallerCallback(typename EnergySignaller::Client* client,
                                                                             Args&&... args)
{
    return client->registerEnergyCallback(std::forward<Args>(args)...);
}

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_SIGNALLERS_H
