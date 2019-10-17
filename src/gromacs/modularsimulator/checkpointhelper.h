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
 * \brief Declares the checkpoint helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_CHECKPOINTHELPER_H
#define GMX_MODULARSIMULATOR_CHECKPOINTHELPER_H

#include <vector>

#include "gromacs/mdlib/checkpointhandler.h"

#include "modularsimulatorinterfaces.h"

struct gmx_walltime_accounting;
struct ObservablesHistory;

namespace gmx
{
class MDLogger;
class TrajectoryElement;

/*! \libinternal
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
 * NS step, the checkpoint helper on master rank signals to all other ranks
 * that checkpointing is about to occur. At the next NS step, the checkpoint
 * is written. On the last step, checkpointing happens immediately after the
 * step (no signalling). To be able to react to last step being signalled,
 * the CheckpointHelper does also implement the `ISimulatorElement` interface,
 * but does only register a function if the last step has been called. It
 * should be placed on top of the simulator loop.
 *
 * Checkpointing happens at the end of a simulation step, which gives a
 * straightforward re-entry point at the top of the simulator loop.
 *
 * In the current implementation, the clients of CheckpointHelper fill a
 * legacy t_state object (passed via pointer) with whatever data they need
 * to store. The CheckpointHelper then writes the t_state object to file.
 * This is an intermediate state of the code, as the long-term plan is for
 * modules to read and write from a checkpoint file directly, without the
 * need for a central object. The current implementation allows, however,
 * to define clearly which modules take part in checkpointing, while using
 * the current infrastructure for reading and writing to checkpoint.
 *
 * \todo Develop this into a module solely providing a file handler to
 *       modules for checkpoint reading and writing.
 */
class CheckpointHelper final : public ILastStepSignallerClient, public ISimulatorElement
{
public:
    //! Constructor
    CheckpointHelper(std::vector<ICheckpointHelperClient*> clients,
                     std::unique_ptr<CheckpointHandler>    checkpointHandler,
                     int                                   initStep,
                     TrajectoryElement*                    trajectoryElement,
                     int                                   globalNumAtoms,
                     FILE*                                 fplog,
                     t_commrec*                            cr,
                     ObservablesHistory*                   observablesHistory,
                     gmx_walltime_accounting*              walltime_accounting,
                     t_state*                              state_global,
                     bool                                  writeFinalCheckpoint);

    /*! \brief Run checkpointing
     *
     * Sets signal and / or performs checkpointing at neighbor searching steps
     *
     * @param step  The step number
     * @param time  The time
     */
    void run(Step step, Time time);

    /*! \brief Register run function for step / time
     *
     * Performs checkpointing at the last step. This is part of the element call
     * list, as the checkpoint helper need to be able to react to the last step
     * being signalled.
     *
     * @param step                 The step number
     * @param time                 The time
     * @param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction) override;

    //! No element setup needed
    void elementSetup() override {}
    //! No element teardown needed
    void elementTeardown() override {}

private:
    //! List of checkpoint clients
    std::vector<ICheckpointHelperClient*> clients_;

    //! The checkpoint handler
    std::unique_ptr<CheckpointHandler> checkpointHandler_;

    //! The first step of the simulation
    const Step initStep_;
    //! The last step of the simulation
    Step lastStep_;
    //! The total number of atoms
    const int globalNumAtoms_;
    //! Whether a checkpoint is written on the last step
    const bool writeFinalCheckpoint_;

    //! ILastStepSignallerClient implementation
    SignallerCallbackPtr registerLastStepCallback() override;

    //! The actual checkpoint writing function
    void writeCheckpoint(Step step, Time time);

    //! Pointer to the trajectory element - to use file pointer
    TrajectoryElement* trajectoryElement_;

    //! A local t_state object to gather data in
    //! {
    std::unique_ptr<t_state> localState_;
    t_state*                 localStateInstance_;
    //! }

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    t_commrec* cr_;
    //! History of simulation observables.
    ObservablesHistory* observablesHistory_;
    //! Manages wall time accounting.
    gmx_walltime_accounting* walltime_accounting_;
    //! Full simulation state (only non-nullptr on master rank).
    t_state* state_global_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_CHECKPOINTHELPER_H
