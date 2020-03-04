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
 * \brief Declares the trajectory element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_TRAJECTORYELEMENT_H
#define GMX_MODULARSIMULATOR_TRAJECTORYELEMENT_H

#include <vector>

#include "gromacs/compat/pointers.h"

#include "modularsimulatorinterfaces.h"

struct gmx_mtop_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_filenm;
struct t_inputrec;

namespace gmx
{
class IMDOutputProvider;
struct MdModulesNotifier;
struct MdrunOptions;
enum class StartingBehavior;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Trajectory element signals and handles trajectory writing
 *
 * The trajectory element is both a signaller and a simulator element.
 *
 * During signalling phase, it checks whether the current step is a writing
 * step for either the energy or the state (position, velocity, forces)
 * trajectory. It then notifies the signaller clients of the upcoming step.
 *
 * For the simulator run, the element registers a run function at trajectory
 * writing steps. Trajectory writing is done using a client system - the
 * element only prepares the output struct, and passes it to the clients who
 * write their part of the trajectory.
 */
class TrajectoryElement final :
    public ISimulatorElement,
    public ISignaller,
    public ILastStepSignallerClient,
    public ILoggingSignallerClient
{
public:
    friend class TrajectoryElementBuilder;

    /*
     * Methods for the signaller part of the element
     */

    /*! \brief Prepare signaller
     *
     * Check that necessary registration was done
     */
    void signallerSetup() override;

    /*! \brief Run the signaller at a specific step / time
     *
     * Informs clients when energy or state will be written.
     *
     * @param step           The current time step
     * @param time           The current time
     */
    void signal(Step step, Time time) override;

    /*
     * Methods for the trajectory writing part of the element
     */

    /*! \brief Prepare trajectory writer
     *
     * During setup, the trajectory writer will query the writer clients for
     * their callbacks. It will also call the setup methods of the different
     * clients. To be run before the main simulator run, but after all clients
     * were registered.
     */
    void elementSetup() override;

    /*! \brief Register run function for step / time
     *
     * Registers a trajectory writing function if the current step / time is
     * either a state or energy writing step, as defined by the signaller
     *
     * @param step                 The step number
     * @param time                 The time
     * @param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction) override;

    /*! \brief Teardown trajectory writer
     *
     * During teardown, the trajectory writer will call the teardown
     * methods of the clients and perform some additional clean-up.
     * To be run after the main simulator run.
     */
    void elementTeardown() override;

    //! @cond
    // (doxygen doesn't like these...)
    //! Allow CheckpointHelper to use outf_ (TODO: Can we improve this?)
    friend class CheckpointHelper;
    //! @endcond

private:
    //! Constructor
    TrajectoryElement(std::vector<SignallerCallbackPtr>     signalEnergyCallbacks,
                      std::vector<SignallerCallbackPtr>     signalStateCallbacks,
                      std::vector<ITrajectoryWriterClient*> writerClients,
                      FILE*                                 fplog,
                      int                                   nfile,
                      const t_filenm                        fnm[],
                      const MdrunOptions&                   mdrunOptions,
                      const t_commrec*                      cr,
                      IMDOutputProvider*                    outputProvider,
                      const MdModulesNotifier&              mdModulesNotifier,
                      const t_inputrec*                     inputrec,
                      gmx_mtop_t*                           top_global,
                      const gmx_output_env_t*               oenv,
                      gmx_wallcycle*                        wcycle,
                      StartingBehavior                      startingBehavior,
                      bool                                  simulationsSharingState);

    //! The next energy writing step
    Step writeEnergyStep_;
    //! The next state writing step
    Step writeStateStep_;
    //! The next communicated log writing step
    Step logWritingStep_;

    //! The output object
    gmx_mdoutf* outf_;

    //! ILoggingSignallerClient implementation
    SignallerCallbackPtr registerLoggingCallback() override;

    /*
     * Signaller
     */
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
    std::vector<SignallerCallbackPtr> signalEnergyCallbacks_;
    std::vector<SignallerCallbackPtr> signalStateCallbacks_;
    //! }

    /*
     * Last step client
     */
    Step lastStep_;
    bool lastStepRegistrationDone_;
    //! ILastStepSignallerClient implementation
    SignallerCallbackPtr registerLastStepCallback() override;

    /*
     * Trajectory writing
     */
    //! The trajectory writing clients
    std::vector<ITrajectoryWriterClient*> writerClients_;

    //! Callbacks to write trajectory
    //! {
    std::vector<ITrajectoryWriterCallbackPtr> runStateCallbacks_;
    std::vector<ITrajectoryWriterCallbackPtr> runEnergyCallbacks_;
    //! }

    //! The writing function - calls the clients to get their contributions
    void write(Step step, Time time, bool writeState, bool writeEnergy, bool writeLog);
};

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Build the `TrajectoryElement`
 *
 * This builder allows clients to register with the trajectory element, either
 * as signaller clients or as writer clients. The builder then builds the
 * element.
 */
class TrajectoryElementBuilder final
{
public:
    //! Allows clients to register to the signaller
    void registerSignallerClient(compat::not_null<ITrajectorySignallerClient*> client);

    //! Allows clients to register as trajectory writers
    void registerWriterClient(compat::not_null<ITrajectoryWriterClient*> client);

    //! Build the TrajectoryElement
    template<typename... Args>
    std::unique_ptr<TrajectoryElement> build(Args&&... args);

private:
    //! List of signaller clients
    std::vector<ITrajectorySignallerClient*> signallerClients_;
    //! List of writer clients
    std::vector<ITrajectoryWriterClient*> writerClients_;
};

template<typename... Args>
std::unique_ptr<TrajectoryElement> TrajectoryElementBuilder::build(Args&&... args)
{
    std::vector<SignallerCallbackPtr> signalEnergyCallbacks;
    std::vector<SignallerCallbackPtr> signalStateCallbacks;
    // Allow clients to register their callbacks
    for (auto& client : signallerClients_)
    {
        // don't register nullptr
        if (auto energyCallback =
                    client->registerTrajectorySignallerCallback(TrajectoryEvent::EnergyWritingStep))
        {
            signalEnergyCallbacks.emplace_back(std::move(energyCallback));
        }
        if (auto stateCallback =
                    client->registerTrajectorySignallerCallback(TrajectoryEvent::StateWritingStep))
        {
            signalStateCallbacks.emplace_back(std::move(stateCallback));
        }
    }
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<TrajectoryElement>(
            new TrajectoryElement(std::move(signalEnergyCallbacks), std::move(signalStateCallbacks),
                                  std::move(writerClients_), std::forward<Args>(args)...));
}

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_TRAJECTORYELEMENT_H
