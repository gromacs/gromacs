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
/*! \internal \file
 * \brief Defines the trajectory element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */
#include "gmxpre.h"

#include "trajectoryelement.h"

#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/inputrec.h"

namespace gmx
{
TrajectoryElement::TrajectoryElement(std::vector<SignallerCallbackPtr>     signalEnergyCallbacks,
                                     std::vector<SignallerCallbackPtr>     signalStateCallbacks,
                                     std::vector<ITrajectoryWriterClient*> writerClients,
                                     FILE*                                 fplog,
                                     int                                   nfile,
                                     const t_filenm                        fnm[],
                                     const MdrunOptions&                   mdrunOptions,
                                     const t_commrec*                      cr,
                                     gmx::IMDOutputProvider*               outputProvider,
                                     const MdModulesNotifier&              mdModulesNotifier,
                                     const t_inputrec*                     inputrec,
                                     gmx_mtop_t*                           top_global,
                                     const gmx_output_env_t*               oenv,
                                     gmx_wallcycle*                        wcycle,
                                     StartingBehavior                      startingBehavior,
                                     const bool                            simulationsShareState) :
    writeEnergyStep_(-1),
    writeStateStep_(-1),
    outf_(init_mdoutf(fplog,
                      nfile,
                      fnm,
                      mdrunOptions,
                      cr,
                      outputProvider,
                      mdModulesNotifier,
                      inputrec,
                      top_global,
                      oenv,
                      wcycle,
                      startingBehavior,
                      simulationsShareState,
                      nullptr)),
    nstxout_(inputrec->nstxout),
    nstvout_(inputrec->nstvout),
    nstfout_(inputrec->nstfout),
    nstxoutCompressed_(inputrec->nstxout_compressed),
    tngBoxOut_(mdoutf_get_tng_box_output_interval(outf_)),
    tngLambdaOut_(mdoutf_get_tng_lambda_output_interval(outf_)),
    tngBoxOutCompressed_(mdoutf_get_tng_compressed_box_output_interval(outf_)),
    tngLambdaOutCompressed_(mdoutf_get_tng_compressed_lambda_output_interval(outf_)),
    nstenergy_(inputrec->nstenergy),
    signalEnergyCallbacks_(std::move(signalEnergyCallbacks)),
    signalStateCallbacks_(std::move(signalStateCallbacks)),
    lastStep_(-1),
    lastStepRegistrationDone_(false),
    writerClients_(std::move(writerClients))
{
}

void TrajectoryElement::signallerSetup()
{
    GMX_ASSERT(lastStepRegistrationDone_,
               "TrajectoryElement needs to be registered to LastStepSignaller.");
}

void TrajectoryElement::signal(Step step, Time time)
{
    if (do_per_step(step, nstxout_) || do_per_step(step, nstvout_) || do_per_step(step, nstfout_)
        || do_per_step(step, nstxoutCompressed_) || do_per_step(step, tngBoxOut_)
        || do_per_step(step, tngLambdaOut_) || do_per_step(step, tngBoxOutCompressed_)
        || do_per_step(step, tngLambdaOutCompressed_))
    {
        writeStateStep_ = step;
        for (const auto& callback : signalStateCallbacks_)
        {
            (*callback)(step, time);
        }
    }

    if (do_per_step(step, nstenergy_) || step == lastStep_)
    {
        writeEnergyStep_ = step;
        for (const auto& callback : signalEnergyCallbacks_)
        {
            (*callback)(step, time);
        }
    }
}

void TrajectoryElement::elementSetup()
{
    for (auto& client : writerClients_)
    {
        auto callback = client->registerTrajectoryWriterCallback(TrajectoryEvent::StateWritingStep);
        if (callback)
        {
            runStateCallbacks_.emplace_back(std::move(callback));
        }
        callback = client->registerTrajectoryWriterCallback(TrajectoryEvent::EnergyWritingStep);
        if (callback)
        {
            runEnergyCallbacks_.emplace_back(std::move(callback));
        }
        client->trajectoryWriterSetup(outf_);
    }
}

void TrajectoryElement::scheduleTask(Step step, Time time, const RegisterRunFunctionPtr& registerRunFunction)
{
    const bool writeEnergyThisStep = writeEnergyStep_ == step;
    const bool writeStateThisStep  = writeStateStep_ == step;
    const bool writeLogThisStep    = logWritingStep_ == step;
    if (writeEnergyThisStep || writeStateThisStep || writeLogThisStep)
    {
        (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
                [this, step, time, writeStateThisStep, writeEnergyThisStep, writeLogThisStep]() {
                    write(step, time, writeStateThisStep, writeEnergyThisStep, writeLogThisStep);
                }));
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
            (*callback)(outf_, step, time, writeState, writeLog);
        }
    }
    if (writeEnergy || writeLog)
    {
        for (auto& callback : runEnergyCallbacks_)
        {
            (*callback)(outf_, step, time, writeEnergy, writeLog);
        }
    }
}

SignallerCallbackPtr TrajectoryElement::registerLastStepCallback()
{
    lastStepRegistrationDone_ = true;
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->lastStep_ = step; });
}

SignallerCallbackPtr TrajectoryElement::registerLoggingCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time /*unused*/) { logWritingStep_ = step; });
}

void TrajectoryElementBuilder::registerSignallerClient(compat::not_null<ITrajectorySignallerClient*> client)
{
    signallerClients_.emplace_back(client);
}

void TrajectoryElementBuilder::registerWriterClient(compat::not_null<ITrajectoryWriterClient*> client)
{
    writerClients_.emplace_back(client);
}

} // namespace gmx
