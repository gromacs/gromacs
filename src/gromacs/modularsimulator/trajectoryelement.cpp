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
 * \brief Defines the trajectory element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */
#include "gmxpre.h"

#include "trajectoryelement.h"

#include <functional>

#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"

namespace gmx
{
TrajectoryElement::TrajectoryElement(std::vector<ITrajectoryWriterClient*> writerClients,
                                     FILE*                                 fplog,
                                     int                                   nfile,
                                     const t_filenm                        fnm[],
                                     const MdrunOptions&                   mdrunOptions,
                                     const t_commrec*                      cr,
                                     gmx::IMDOutputProvider*               outputProvider,
                                     const MDModulesNotifiers&             mdModulesNotifiers,
                                     const t_inputrec*                     inputrec,
                                     const gmx_mtop_t&                     top_global,
                                     const gmx_output_env_t*               oenv,
                                     gmx_wallcycle*                        wcycle,
                                     StartingBehavior                      startingBehavior,
                                     const bool                            simulationsShareState) :
    writeEnergyStep_(-1),
    writeStateStep_(-1),
    writeLogStep_(-1),
    outf_(init_mdoutf(fplog,
                      nfile,
                      fnm,
                      mdrunOptions,
                      cr,
                      outputProvider,
                      mdModulesNotifiers,
                      inputrec,
                      top_global,
                      oenv,
                      wcycle,
                      startingBehavior,
                      simulationsShareState,
                      nullptr)),
    writerClients_(std::move(writerClients))
{
}

int TrajectoryElement::tngBoxOut() const
{
    return mdoutf_get_tng_box_output_interval(outf_);
}
int TrajectoryElement::tngLambdaOut() const
{
    return mdoutf_get_tng_lambda_output_interval(outf_);
}
int TrajectoryElement::tngBoxOutCompressed() const
{
    return mdoutf_get_tng_compressed_box_output_interval(outf_);
}
int TrajectoryElement::tngLambdaOutCompressed() const
{
    return mdoutf_get_tng_compressed_lambda_output_interval(outf_);
}

void TrajectoryElement::elementSetup()
{
    for (auto& client : writerClients_)
    {
        auto callback = client->registerTrajectoryWriterCallback(TrajectoryEvent::StateWritingStep);
        if (callback)
        {
            runStateCallbacks_.emplace_back(std::move(*callback));
        }
        callback = client->registerTrajectoryWriterCallback(TrajectoryEvent::EnergyWritingStep);
        if (callback)
        {
            runEnergyCallbacks_.emplace_back(std::move(*callback));
        }
        client->trajectoryWriterSetup(outf_);
    }
}

void TrajectoryElement::scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction)
{
    const bool writeEnergyThisStep = writeEnergyStep_ == step;
    const bool writeStateThisStep  = writeStateStep_ == step;
    const bool writeLogThisStep    = writeLogStep_ == step;
    if (writeEnergyThisStep || writeStateThisStep || writeLogThisStep)
    {
        registerRunFunction([this, step, time, writeStateThisStep, writeEnergyThisStep, writeLogThisStep]() {
            write(step, time, writeStateThisStep, writeEnergyThisStep, writeLogThisStep);
        });
    }
}

void TrajectoryElement::elementTeardown()
{
    for (auto& client : writerClients_)
    {
        client->trajectoryWriterTeardown(outf_);
    }
    mdoutf_tng_close(outf_);
    done_mdoutf(outf_);
}

void TrajectoryElement::write(Step step, Time time, bool writeState, bool writeEnergy, bool writeLog)
{
    if (writeState || writeLog)
    {
        for (auto& callback : runStateCallbacks_)
        {
            callback(outf_, step, time, writeState, writeLog);
        }
    }
    if (writeEnergy || writeLog)
    {
        for (auto& callback : runEnergyCallbacks_)
        {
            callback(outf_, step, time, writeEnergy, writeLog);
        }
    }
}

std::optional<SignallerCallback> TrajectoryElement::registerLoggingCallback()
{
    return [this](Step step, Time /*unused*/) { writeLogStep_ = step; };
}

std::optional<SignallerCallback> TrajectoryElement::registerTrajectorySignallerCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::StateWritingStep)
    {
        return [this](Step step, Time /*unused*/) { this->writeStateStep_ = step; };
    }
    if (event == TrajectoryEvent::EnergyWritingStep)
    {
        return [this](Step step, Time /*unused*/) { this->writeEnergyStep_ = step; };
    }
    return std::nullopt;
}

void TrajectoryElementBuilder::registerWriterClient(ITrajectoryWriterClient* client)
{
    if (client)
    {
        if (state_ == ModularSimulatorBuilderState::NotAcceptingClientRegistrations)
        {
            throw SimulationAlgorithmSetupError(
                    "Tried to register to signaller after it was built.");
        }
        writerClients_.emplace_back(client);
    }
}

} // namespace gmx
