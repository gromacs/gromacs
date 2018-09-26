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

#include <cmath>

#include "gromacs/timing/walltime_accounting.h"

namespace gmx
{

CheckpointHandler::CheckpointHandler(
        compat::not_null<AccumulateGlobalsBuilder*> accumulateGlobalsBuilder,
        bool                                        simulationsShareState,
        bool                                        neverUpdateNeighborList,
        bool                                        isMaster,
        bool                                        writeFinalCheckpoint,
        real                                        checkpointingPeriod,
        int                                         numMpiThreads) :
    hasUnhandledSignal_(false),
    waitingForSignalReduction_(false),
    checkpointThisStep_(false),
    numberOfNextCheckpoint_(1),
    rankCanSetSignal_(checkpointingPeriod >= 0 && isMaster),
    checkpointingIsActive_(checkpointingPeriod >= 0),
    writeFinalCheckpoint_(writeFinalCheckpoint),
    neverUpdateNeighborlist_(neverUpdateNeighborList),
    checkpointingPeriod_(checkpointingPeriod),
    numMpiThreads_(numMpiThreads)
{
    if (checkpointingIsActive_)
    {
        accumulateGlobalsBuilder->registerClient(
                compat::not_null<IAccumulateGlobalsClient *>(this), simulationsShareState);
    }
}

void CheckpointHandler::setSignalImpl(
        gmx_walltime_accounting_t walltime_accounting) const
{
    if (!hasUnhandledSignal_ && !waitingForSignalReduction_ &&
        (checkpointingPeriod_ == 0 ||
         walltime_accounting_get_time_since_start(walltime_accounting) >=
         numberOfNextCheckpoint_ * checkpointingPeriod_ * 60.0))
    {
        signalView_[0]             = static_cast<double>(CheckpointSignal::doCheckpoint);
        waitingForSignalReduction_ = true;
    }
}

void CheckpointHandler::decideIfCheckpointingThisStepImpl(
        bool bNS, bool bFirstStep, bool bLastStep)
{
    if (static_cast<CheckpointSignal>(lround(signalView_[0])) == CheckpointSignal::doCheckpoint &&
        signalView_[1] == numMpiThreads_)
    {
        signalView_[0]             = static_cast<double>(CheckpointSignal::noSignal);
        hasUnhandledSignal_        = true;
        waitingForSignalReduction_ = false;
    }
    signalView_[1] = 1;

    if (((hasUnhandledSignal_ && (bNS || neverUpdateNeighborlist_)) ||
         (bLastStep && writeFinalCheckpoint_)) &&
        !bFirstStep)
    {
        checkpointThisStep_ = true;
        hasUnhandledSignal_ = false;
        numberOfNextCheckpoint_++;
    }
}

int CheckpointHandler::getNumGlobalsRequired() const
{
    return 2;
}

void CheckpointHandler::setViewForGlobals(
        AccumulateGlobals *accumulateGlobals gmx_unused,
        ArrayRef<double>   view)
{
    // saving the ref to accumulateGlobals would only be interesting if we
    // needed to signal that global reduction is needed - but these signals
    // are not urgent.

    signalView_ = view;
    // set signal empty
    signalView_[0] = static_cast<double>(CheckpointSignal::noSignal);
    // The second reduced variable is notifying that reduction has happened -
    // we set it to zero on all ranks when we set the signal, and it will
    // hence be reduced to the number of ranks once reduction has happened.
    signalView_[1] = checkpointingIsActive_ ? 1 : 0;
}

void CheckpointHandler::notifyAfterCommunication()
{
    // Not sure if we'll need that for anything...
}

} // namespace gmx
