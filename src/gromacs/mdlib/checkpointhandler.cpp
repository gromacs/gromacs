/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief
 * Defines the checkpoint handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "checkpointhandler.h"

#include "gromacs/timing/walltime_accounting.h"

namespace gmx
{

CheckpointHandler::CheckpointHandler(
        compat::not_null<AccumulatorBuilder<ISimulationAccumulatorClient>*>      simulationAccumulatorBuilder,
        compat::not_null<AccumulatorBuilder<IMultiSimulationAccumulatorClient>*> multiSimulationAccumulatorBuilder,
        bool                                                                     neverUpdateNeighborList,
        bool                                                                     isMaster,
        bool                                                                     writeFinalCheckpoint,
        real                                                                     checkpointingPeriod) :
    multiSimAccumulator_(nullptr),
    doMultiSim_(false),
    hasUnhandledSignal_(false),
    checkpointThisStep_(false),
    numberOfNextCheckpoint_(1),
    rankCanSetSignal_(checkpointingPeriod >= 0 && isMaster),
    checkpointingIsActive_(checkpointingPeriod >= 0),
    writeFinalCheckpoint_(writeFinalCheckpoint),
    neverUpdateNeighborlist_(neverUpdateNeighborList),
    checkpointingPeriod_(checkpointingPeriod)
{
    if (checkpointingIsActive_)
    {
        simulationAccumulatorBuilder->registerClient(
                compat::not_null<ISimulationAccumulatorClient*>(this));
        multiSimulationAccumulatorBuilder->registerClient(
                compat::not_null<IMultiSimulationAccumulatorClient*>(this));

    }
}

void CheckpointHandler::setSignalImpl(
        gmx_walltime_accounting_t walltime_accounting)
{
    const bool signalNeeded = !hasUnhandledSignal_ &&
        (checkpointingPeriod_ == 0 ||
         walltime_accounting_get_time_since_start(walltime_accounting) >=
         numberOfNextCheckpoint_ * checkpointingPeriod_ * 60.0);
    const bool multiSimSignalSet =
        doMultiSim_ && toCheckpointSignal(multiSimSignalView_[0]) == CheckpointSignal::doCheckpoint;
    if (signalNeeded && doMultiSim_ && !multiSimSignalSet)
    {
        multiSimSignalView_[0] = static_cast<double>(CheckpointSignal::doCheckpoint);
        multiSimAccumulator_->notifyReductionRequired(compat::not_null<IMultiSimulationAccumulatorClient*>(this));
    }
    else if (signalNeeded || multiSimSignalSet)
    {
        signalView_[0] = static_cast<double>(CheckpointSignal::doCheckpoint);
    }
}

void CheckpointHandler::decideIfCheckpointingThisStepImpl(
        bool bNS, bool bFirstStep, bool bLastStep)
{
    if (((hasUnhandledSignal_ && (bNS || neverUpdateNeighborlist_)) ||
         (bLastStep && writeFinalCheckpoint_)) &&
        !bFirstStep)
    {
        checkpointThisStep_ = true;
        hasUnhandledSignal_ = false;
        numberOfNextCheckpoint_++;
    }
}

int CheckpointHandler::getNumSimulationGlobalsRequired() const
{
    return 1;
}

void CheckpointHandler::setViewForSimulationGlobals(
        Accumulator<ISimulationAccumulatorClient> *accumulator gmx_unused,
        ArrayRef<double>                           view)
{
    // No need for SimulationAccumulator pointer

    signalView_ = view;
    // set signal empty
    signalView_[0] = static_cast<double>(CheckpointSignal::noSignal);
}

void CheckpointHandler::notifyAfterSimulationCommunication()
{
    // Not sure if we'll need that for anything...
}

int CheckpointHandler::getNumMultiSimulationGlobalsRequired() const
{
    return 1;
}

void CheckpointHandler::setViewForMultiSimulationGlobals(
        Accumulator<IMultiSimulationAccumulatorClient> *accumulator,
        ArrayRef<double>                                view)
{
    multiSimAccumulator_ = accumulator;
    doMultiSim_          = true;
    multiSimSignalView_  = view;
    // set signal empty
    multiSimSignalView_[0] = static_cast<double>(CheckpointSignal::noSignal);
}

void CheckpointHandler::notifyAfterMultiSimulationCommunication()
{
    // Not sure if we'll need that for anything...
}

} // namespace gmx
