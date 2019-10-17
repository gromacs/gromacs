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
/*! \internal \file
 * \brief Defines the checkpoint helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "checkpointhelper.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/state.h"

#include "trajectoryelement.h"

namespace gmx
{
CheckpointHelper::CheckpointHelper(std::vector<ICheckpointHelperClient*> clients,
                                   std::unique_ptr<CheckpointHandler>    checkpointHandler,
                                   int                                   initStep,
                                   TrajectoryElement*                    trajectoryElement,
                                   int                                   globalNumAtoms,
                                   FILE*                                 fplog,
                                   t_commrec*                            cr,
                                   ObservablesHistory*                   observablesHistory,
                                   gmx_walltime_accounting*              walltime_accounting,
                                   t_state*                              state_global,
                                   bool                                  writeFinalCheckpoint) :
    clients_(std::move(clients)),
    checkpointHandler_(std::move(checkpointHandler)),
    initStep_(initStep),
    lastStep_(-1),
    globalNumAtoms_(globalNumAtoms),
    writeFinalCheckpoint_(writeFinalCheckpoint),
    trajectoryElement_(trajectoryElement),
    localState_(nullptr),
    fplog_(fplog),
    cr_(cr),
    observablesHistory_(observablesHistory),
    walltime_accounting_(walltime_accounting),
    state_global_(state_global)
{
    // Get rid of nullptr in clients list
    clients_.erase(std::remove_if(clients_.begin(), clients_.end(),
                                  [](ICheckpointHelperClient* ptr) { return ptr == nullptr; }),
                   clients_.end());
    if (DOMAINDECOMP(cr))
    {
        localState_ = std::make_unique<t_state>();
        dd_init_local_state(cr->dd, state_global, localState_.get());
        localStateInstance_ = localState_.get();
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        localStateInstance_ = state_global;
    }
}

void CheckpointHelper::run(Step step, Time time)
{
    // reads out signal, decides if we should signal checkpoint
    checkpointHandler_->decideIfCheckpointingThisStep(true, step == initStep_, false);
    if (checkpointHandler_->isCheckpointingStep())
    {
        writeCheckpoint(step, time);
    }

    // decides if we should set a checkpointing signal
    checkpointHandler_->setSignal(walltime_accounting_);
}

void CheckpointHelper::scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction)
{
    // Only last step checkpointing is done here
    if (step != lastStep_ || !writeFinalCheckpoint_)
    {
        return;
    }
    (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
            [this, step, time]() { writeCheckpoint(step, time); }));
}

void CheckpointHelper::writeCheckpoint(Step step, Time time)
{
    localStateInstance_->flags = 0;
    for (const auto& client : clients_)
    {
        client->writeCheckpoint(localStateInstance_, state_global_);
    }

    mdoutf_write_to_trajectory_files(fplog_, cr_, trajectoryElement_->outf_, MDOF_CPT,
                                     globalNumAtoms_, step, time, localStateInstance_,
                                     state_global_, observablesHistory_, ArrayRef<RVec>());
}

SignallerCallbackPtr CheckpointHelper::registerLastStepCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->lastStep_ = step; });
}

} // namespace gmx
