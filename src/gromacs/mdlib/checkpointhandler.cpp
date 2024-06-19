/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief
 * Defines the checkpoint handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "checkpointhandler.h"

#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Convert signed char (as used by SimulationSignal) to CheckpointSignal enum
 *
 * Expected values are
 *   \p sig == 0 -- no signal
 *   \p sig >= 1 -- signal received
 */
static inline CheckpointSignal convertToCheckpointSignal(signed char sig)
{
    GMX_ASSERT(sig >= 0, "Unexpected checkpoint signal < 0 received");
    return sig >= 1 ? CheckpointSignal::doCheckpoint : CheckpointSignal::noSignal;
}

CheckpointHandler::CheckpointHandler(compat::not_null<SimulationSignal*> signal,
                                     bool                                simulationsShareState,
                                     bool                                neverUpdateNeighborList,
                                     bool                                isMain,
                                     bool                                writeFinalCheckpoint,
                                     real                                checkpointingPeriod) :
    signal_(*signal),
    checkpointThisStep_(false),
    numberOfNextCheckpoint_(1),
    rankCanSetSignal_(checkpointingPeriod >= 0 && isMain),
    checkpointingIsActive_(checkpointingPeriod >= 0),
    writeFinalCheckpoint_(writeFinalCheckpoint),
    neverUpdateNeighborlist_(neverUpdateNeighborList),
    checkpointingPeriod_(checkpointingPeriod)
{
    if (simulationsShareState)
    {
        signal_.isLocal = false;
    }
}

void CheckpointHandler::setSignalImpl(gmx_walltime_accounting_t walltime_accounting) const
{
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if (convertToCheckpointSignal(signal_.set) == CheckpointSignal::noSignal
        && convertToCheckpointSignal(signal_.sig) == CheckpointSignal::noSignal
        && (checkpointingPeriod_ == 0
            || secondsSinceStart >= numberOfNextCheckpoint_ * checkpointingPeriod_ * 60.0))
    {
        signal_.sig = static_cast<signed char>(CheckpointSignal::doCheckpoint);
    }
}

void CheckpointHandler::decideIfCheckpointingThisStepImpl(bool bNS, bool bFirstStep, bool bLastStep)
{
    checkpointThisStep_ = (((convertToCheckpointSignal(signal_.set) == CheckpointSignal::doCheckpoint
                             && (bNS || neverUpdateNeighborlist_))
                            || (bLastStep && writeFinalCheckpoint_))
                           && !bFirstStep);
    if (checkpointThisStep_)
    {
        signal_.set = static_cast<signed char>(CheckpointSignal::noSignal);
        numberOfNextCheckpoint_++;
    }
}

} // namespace gmx
