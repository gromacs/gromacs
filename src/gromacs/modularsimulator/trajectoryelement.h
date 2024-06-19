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
 * \brief Declares the trajectory element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_TRAJECTORYELEMENT_H
#define GMX_MODULARSIMULATOR_TRAJECTORYELEMENT_H

#include <cstdio>

#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "gromacs/compat/pointers.h"

#include "modularsimulatorinterfaces.h"

struct gmx_mtop_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct gmx_mdoutf;

namespace gmx
{
class IMDOutputProvider;
struct MDModulesNotifiers;
struct MdrunOptions;
enum class StartingBehavior;

/*! \internal
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
class TrajectoryElement final : public ISimulatorElement, public ILoggingSignallerClient, public ITrajectorySignallerClient
{
public:
    friend class TrajectoryElementBuilder;
    //! Get the box writeout frequency for TNG
    [[nodiscard]] int tngBoxOut() const;
    //! Get the lambda writeout frequency for TNG
    [[nodiscard]] int tngLambdaOut() const;
    //! Get the compressed box writeout frequency for TNG
    [[nodiscard]] int tngBoxOutCompressed() const;
    //! Get the compressed lambda writeout frequency for TNG
    [[nodiscard]] int tngLambdaOutCompressed() const;

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
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    /*! \brief Teardown trajectory writer
     *
     * During teardown, the trajectory writer will call the teardown
     * methods of the clients and perform some additional clean-up.
     * To be run after the main simulator run.
     */
    void elementTeardown() override;

    //! \cond
    // (doxygen doesn't like these...)
    //! Allow CheckpointHelper to use outf_ (TODO: Can we improve this?)
    friend class CheckpointHelper;
    //! \endcond

private:
    //! Constructor
    TrajectoryElement(std::vector<ITrajectoryWriterClient*> writerClients,
                      FILE*                                 fplog,
                      int                                   nfile,
                      const t_filenm                        fnm[],
                      const MdrunOptions&                   mdrunOptions,
                      const t_commrec*                      cr,
                      IMDOutputProvider*                    outputProvider,
                      const MDModulesNotifiers&             mdModulesNotifiers,
                      const t_inputrec*                     inputrec,
                      const gmx_mtop_t&                     top_global,
                      const gmx_output_env_t*               oenv,
                      gmx_wallcycle*                        wcycle,
                      StartingBehavior                      startingBehavior,
                      bool                                  simulationsSharingState);

    //! The next energy writing step
    Step writeEnergyStep_;
    //! The next state writing step
    Step writeStateStep_;
    //! The next communicated log writing step
    Step writeLogStep_;

    //! The output object
    gmx_mdoutf* outf_;

    //! ILoggingSignallerClient implementation
    std::optional<SignallerCallback> registerLoggingCallback() override;
    //! ITrajectorySignallerClient implementation
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;

    /*
     * Trajectory writing
     */
    //! The trajectory writing clients
    std::vector<ITrajectoryWriterClient*> writerClients_;

    //! Callbacks to write trajectory
    //! {
    std::vector<ITrajectoryWriterCallback> runStateCallbacks_;
    std::vector<ITrajectoryWriterCallback> runEnergyCallbacks_;
    //! }

    //! The writing function - calls the clients to get their contributions
    void write(Step step, Time time, bool writeState, bool writeEnergy, bool writeLog);
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Build the `TrajectoryElement`
 *
 * This builder allows clients to register with the trajectory element
 * as writer clients. The builder then builds the trajectory element.
 */
class TrajectoryElementBuilder final
{
public:
    //! Allows clients to register as trajectory writers
    void registerWriterClient(ITrajectoryWriterClient* client);

    //! Build the TrajectoryElement
    template<typename... Args>
    std::unique_ptr<TrajectoryElement> build(Args&&... args);

private:
    //! List of writer clients
    std::vector<ITrajectoryWriterClient*> writerClients_;
    //! The state of the builder
    ModularSimulatorBuilderState state_ = ModularSimulatorBuilderState::AcceptingClientRegistrations;
};

template<typename... Args>
std::unique_ptr<TrajectoryElement> TrajectoryElementBuilder::build(Args&&... args)
{
    state_ = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<TrajectoryElement>(
            new TrajectoryElement(std::move(writerClients_), std::forward<Args>(args)...));
}

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_TRAJECTORYELEMENT_H
