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
 * \brief Defines the signallers for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */
#include "gmxpre.h"

#include "signallers.h"

#include <algorithm>

#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/stophandler.h"

#include "modularsimulatorinterfaces.h"

namespace gmx
{
//! Helper function to call all callbacks in a list
static inline void runAllCallbacks(std::vector<SignallerCallbackPtr>& callbacks, Step step, Time time)
{
    for (const auto& callback : callbacks)
    {
        (*callback)(step, time);
    }
}

NeighborSearchSignaller::NeighborSearchSignaller(std::vector<SignallerCallbackPtr> callbacks,
                                                 Step                              nstlist,
                                                 Step                              initStep,
                                                 Time                              initTime) :
    callbacks_(std::move(callbacks)),
    nstlist_(nstlist),
    initStep_(initStep),
    initTime_(initTime)
{
}

void NeighborSearchSignaller::signal(Step step, Time time)
{
    // Neighbor search happens at regular intervals, and always on first step of simulation
    if (do_per_step(step, nstlist_) || step == initStep_)
    {
        runAllCallbacks(callbacks_, step, time);
    }
}

LastStepSignaller::LastStepSignaller(std::vector<SignallerCallbackPtr> callbacks,
                                     gmx::Step                         nsteps,
                                     gmx::Step                         initStep,
                                     StopHandler*                      stopHandler) :
    callbacks_(std::move(callbacks)),
    stopStep_(initStep + nsteps),
    signalledStopCondition_(false),
    stopHandler_(stopHandler),
    nextNSStep_(-1),
    nsStepRegistrationDone_(false)
{
}

void LastStepSignaller::signal(Step step, Time time)
{
    if (signalledStopCondition_)
    {
        return;
    }
    bool isNSStep           = (step == nextNSStep_);
    signalledStopCondition_ = stopHandler_->stoppingAfterCurrentStep(isNSStep);
    if (step == stopStep_ || signalledStopCondition_)
    {
        runAllCallbacks(callbacks_, step, time);
    }
}

void LastStepSignaller::signallerSetup()
{
    GMX_ASSERT(nsStepRegistrationDone_,
               "LastStepSignaller needs to be registered to NeighborSearchSignaller.");
}

SignallerCallbackPtr LastStepSignaller::registerNSCallback()
{
    nsStepRegistrationDone_ = true;
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->nextNSStep_ = step; });
}

LoggingSignaller::LoggingSignaller(std::vector<SignallerCallbackPtr> callbacks,
                                   Step                              nstlog,
                                   Step                              initStep,
                                   Time                              initTime) :
    callbacks_(std::move(callbacks)),
    nstlog_(nstlog),
    initStep_(initStep),
    initTime_(initTime),
    lastStep_(-1),
    lastStepRegistrationDone_(false)
{
}

void LoggingSignaller::signal(Step step, Time time)
{
    if (do_per_step(step, nstlog_) || step == lastStep_)
    {
        runAllCallbacks(callbacks_, step, time);
    }
}

void LoggingSignaller::signallerSetup()
{
    GMX_ASSERT(lastStepRegistrationDone_,
               "LoggingSignaller needs to be registered to LastStepSignaller.");
}

SignallerCallbackPtr LoggingSignaller::registerLastStepCallback()
{
    lastStepRegistrationDone_ = true;
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->lastStep_ = step; });
}

EnergySignaller::EnergySignaller(std::vector<SignallerCallbackPtr> calculateEnergyCallbacks,
                                 std::vector<SignallerCallbackPtr> calculateVirialCallbacks,
                                 std::vector<SignallerCallbackPtr> calculateFreeEnergyCallbacks,
                                 int                               nstcalcenergy,
                                 int                               nstcalcfreeenergy,
                                 int                               nstcalcvirial) :
    calculateEnergyCallbacks_(std::move(calculateEnergyCallbacks)),
    calculateVirialCallbacks_(std::move(calculateVirialCallbacks)),
    calculateFreeEnergyCallbacks_(std::move(calculateFreeEnergyCallbacks)),
    nstcalcenergy_(nstcalcenergy),
    nstcalcfreeenergy_(nstcalcfreeenergy),
    nstcalcvirial_(nstcalcvirial),
    energyWritingStep_(-1),
    trajectoryRegistrationDone_(false),
    loggingStep_(-1),
    loggingRegistrationDone_(false)
{
}

void EnergySignaller::signal(Step step, Time time)
{
    bool calculateEnergy     = do_per_step(step, nstcalcenergy_);
    bool calculateFreeEnergy = do_per_step(step, nstcalcfreeenergy_);
    bool calculateVirial     = do_per_step(step, nstcalcvirial_);
    bool writeEnergy         = energyWritingStep_ == step;

    if (calculateEnergy || writeEnergy || step == loggingStep_)
    {
        runAllCallbacks(calculateEnergyCallbacks_, step, time);
    }
    if (calculateEnergy || writeEnergy || step == loggingStep_ || calculateVirial)
    {
        runAllCallbacks(calculateVirialCallbacks_, step, time);
    }
    if (calculateFreeEnergy)
    {
        runAllCallbacks(calculateFreeEnergyCallbacks_, step, time);
    }
}

void EnergySignaller::signallerSetup()
{
    GMX_ASSERT(loggingRegistrationDone_,
               "EnergySignaller needs to be registered to LoggingSignaller.");
    GMX_ASSERT(trajectoryRegistrationDone_,
               "EnergySignaller needs to be registered to TrajectoryElement.");
}

SignallerCallbackPtr EnergySignaller::registerTrajectorySignallerCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::EnergyWritingStep)
    {
        trajectoryRegistrationDone_ = true;
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time gmx_unused time) { this->energyWritingStep_ = step; });
    }
    return nullptr;
}

SignallerCallbackPtr EnergySignaller::registerLoggingCallback()
{
    loggingRegistrationDone_ = true;
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->loggingStep_ = step; });
}

} // namespace gmx
