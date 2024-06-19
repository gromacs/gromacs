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
 * \brief Declares the checkpoint helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_CHECKPOINTHELPER_H
#define GMX_MODULARSIMULATOR_CHECKPOINTHELPER_H

#include <cstdio>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/utility/exceptions.h"

#include "modularsimulatorinterfaces.h"

struct gmx_walltime_accounting;
struct ObservablesHistory;
class t_state;
struct t_commrec;

namespace gmx
{
class KeyValueTreeObject;
class MDLogger;
class TrajectoryElement;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Checkpoint helper
 *
 * The `CheckpointHelper` is responsible to write checkpoints. In the
 * longer term, it will also be responsible to read checkpoints, but this
 * is not yet implemented.
 *
 * Writing checkpoints is done just before neighbor-searching (NS) steps,
 * or after the last step. Checkpointing occurs periodically (by default,
 * every 15 minutes), and needs two NS steps to take effect - on the first
 * NS step, the checkpoint helper on main rank signals to all other ranks
 * that checkpointing is about to occur. At the next NS step, the checkpoint
 * is written. On the last step, checkpointing happens immediately before the
 * step (no signalling). To be able to react to last step being signalled,
 * the CheckpointHelper does also implement the `ISimulatorElement` interface,
 * but does only register a function if the last step has been called. It
 * should be placed on top of the simulator loop.
 *
 * Checkpointing happens at the end of a simulation step, which gives a
 * straightforward re-entry point at the top of the simulator loop.
 *
 * Checkpoint writing is done by passing sub-objects of a
 * WriteCheckpointDataHolder object to the clients. Checkpoint reading is
 * done by passing sub-objects of a ReadCheckpointDataHolder object (passed
 * in from runner level) do the clients.
 *
 * \see ReadCheckpointDataHolder
 * \see WriteCheckpointDataHolder
 * \see CheckpointData
 */
class CheckpointHelper final : public ILastStepSignallerClient, public ISimulatorElement
{
public:
    //! Constructor
    CheckpointHelper(std::vector<std::tuple<std::string, ICheckpointHelperClient*>>&& clients,
                     std::unique_ptr<CheckpointHandler> checkpointHandler,
                     int                                initStep,
                     TrajectoryElement*                 trajectoryElement,
                     FILE*                              fplog,
                     t_commrec*                         cr,
                     ObservablesHistory*                observablesHistory,
                     gmx_walltime_accounting*           walltime_accounting,
                     t_state*                           state_global,
                     bool                               writeFinalCheckpoint);

    /*! \brief Run checkpointing
     *
     * Sets signal and / or performs checkpointing at neighbor searching steps
     *
     * \param step  The step number
     * \param time  The time
     */
    void run(Step step, Time time);

    /*! \brief Register run function for step / time
     *
     * Performs checkpointing at the last step. This is part of the element call
     * list, as the checkpoint helper need to be able to react to the last step
     * being signalled.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! No element setup needed
    void elementSetup() override {}
    //! No element teardown needed
    void elementTeardown() override {}

private:
    //! List of checkpoint clients
    std::vector<std::tuple<std::string, ICheckpointHelperClient*>> clients_;

    //! The checkpoint handler
    std::unique_ptr<CheckpointHandler> checkpointHandler_;

    //! The first step of the simulation
    const Step initStep_;
    //! The last step of the simulation
    Step lastStep_;
    //! Whether a checkpoint is written on the last step
    const bool writeFinalCheckpoint_;

    //! ILastStepSignallerClient implementation
    std::optional<SignallerCallback> registerLastStepCallback() override;

    //! The actual checkpoint writing function
    void writeCheckpoint(Step step, Time time);

    //! Pointer to the trajectory element - to use file pointer
    TrajectoryElement* trajectoryElement_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    t_commrec* cr_;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory_;
    //! Manages wall time accounting.
    gmx_walltime_accounting* walltime_accounting_;
    //! Full simulation state (only non-nullptr on main rank).
    t_state* state_global_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Builder for the checkpoint helper
 */
class CheckpointHelperBuilder
{
public:
    //! Constructor
    CheckpointHelperBuilder(std::unique_ptr<ReadCheckpointDataHolder> checkpointDataHolder,
                            StartingBehavior                          startingBehavior,
                            t_commrec*                                cr);

    //! Register checkpointing client
    void registerClient(ICheckpointHelperClient* client);

    //! Set CheckpointHandler
    void setCheckpointHandler(std::unique_ptr<CheckpointHandler> checkpointHandler);

    //! Return CheckpointHelper
    template<typename... Args>
    std::unique_ptr<CheckpointHelper> build(Args&&... args);

private:
    //! Map of checkpoint clients
    std::map<std::string, ICheckpointHelperClient*> clientsMap_;
    //! Whether we are resetting from checkpoint
    const bool resetFromCheckpoint_;
    //! The input checkpoint data
    std::unique_ptr<ReadCheckpointDataHolder> checkpointDataHolder_;
    //! The checkpoint handler
    std::unique_ptr<CheckpointHandler> checkpointHandler_;
    //! Handles communication.
    t_commrec* cr_;
    //! Whether the builder accepts registrations.
    ModularSimulatorBuilderState state_;
};

template<typename... Args>
std::unique_ptr<CheckpointHelper> CheckpointHelperBuilder::build(Args&&... args)
{
    state_ = ModularSimulatorBuilderState::NotAcceptingClientRegistrations;
    // Make sure that we don't have unused entries in checkpoint
    if (resetFromCheckpoint_)
    {
        for (const auto& key : checkpointDataHolder_->keys())
        {
            if (clientsMap_.count(key) == 0)
            {
                // We have an entry in checkpointDataHolder_ which has no matching client
                throw CheckpointError("Checkpoint entry " + key + " was not read. This "
                                      "likely means that you are not using the same algorithm "
                                      "that was used to create the checkpoint file.");
            }
        }
    }

    std::vector<std::tuple<std::string, ICheckpointHelperClient*>>&& clients = { clientsMap_.begin(),
                                                                                 clientsMap_.end() };
    return std::make_unique<CheckpointHelper>(
            std::move(clients), std::move(checkpointHandler_), std::forward<Args>(args)...);
}

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_CHECKPOINTHELPER_H
