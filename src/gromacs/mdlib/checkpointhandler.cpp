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

#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"

using namespace gmx;

CheckpointHandler::CheckpointHandler(
        compat::not_null<gmx::SimulationSignal*> signal,
        bool                                     simulationsShareState,
        bool                                     neverUpdateNeighborList,
        bool                                     isMaster,
        bool                                     writeFinalCheckpoint,
        real                                     checkpointingPeriod) :
    signal_(*signal),
    rankCanSetSignal_(false),
    checkpointingIsActive_(false),
    checkpointThisStep_(false),
    numberOfNextCheckpoint_(1),
    writeFinalCheckpoint_(writeFinalCheckpoint),
    neverUpdateNeighborlist_(neverUpdateNeighborList),
    checkpointingPeriod_(checkpointingPeriod)
{
    if (simulationsShareState)
    {
        signal_.isLocal = false;
    }

    if (checkpointingPeriod_ >= 0)
    {
        if (isMaster)
        {
            rankCanSetSignal_ = true;
        }
        checkpointingIsActive_ = true;
    }
}

void CheckpointHandler::setSignalImpl_(
        gmx_walltime_accounting_t walltime_accounting)
{
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if (signal_.set == 0 && signal_.sig == 0 &&
        (checkpointingPeriod_ == 0 || secondsSinceStart >= numberOfNextCheckpoint_ * checkpointingPeriod_ * 60.0))
    {
        signal_.sig = 1;
    }
}

void CheckpointHandler::doCheckpointImpl_(
        bool bNS, bool bFirstStep, bool bLastStep)
{
    checkpointThisStep_ = (((signal_.set != 0 && (bNS || neverUpdateNeighborlist_)) ||
                            (bLastStep && writeFinalCheckpoint_)) &&
                           !bFirstStep);
    if (checkpointThisStep_)
    {
        signal_.set = 0;
        numberOfNextCheckpoint_++;
    }
}
