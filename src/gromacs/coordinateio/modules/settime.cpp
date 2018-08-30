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
/*!\file
 * \internal
 * \brief
 * Implements settime class.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "settime.h"

#include <algorithm>

#include "gromacs/utility/exceptions.h"

namespace gmx
{

void
SetTime::processFrame(const int framenumber, t_trxframe *input)
{
    if (!haveProcessedFirstFrame_)
    {
        setInitialTime(input->time);
    }
    if (haveProcessedFirstFrame_ && !haveProcessedSecondFrame_)
    {
        setTimeShift(input->time);
    }

    switch (frameTime_)
    {
        case (ChangeFrameTimeType::efTimeStep):
        case (ChangeFrameTimeType::efBothTime):
            input->time = startTime_ + framenumber*timeShift_;
            break;
        case (ChangeFrameTimeType::efStartTime):
            startTime_ += timeShift_;
            input->time = startTime_;
            break;
        case (ChangeFrameTimeType::efUnchanged):
            break;
        default:
            GMX_THROW(InconsistentInputError("Value for frameTime flag is not supported"));
    }

    input->bTime = true;
}

void
SetTime::setInitialTime(real initialTime)
{
    switch (frameTime_)
    {
        case (ChangeFrameTimeType::efTimeStep):
            startTime_ = initialTime;
            timeShift_ = timeStep_;
            break;
        case (ChangeFrameTimeType::efStartTime):
            // Store old time in timeStep, set initial timeShift_ to 0
            timeStep_  = initialTime;
            timeShift_ = 0;
            break;
        case (ChangeFrameTimeType::efBothTime):
            timeShift_ = timeStep_;
            break;
        case (ChangeFrameTimeType::efUnchanged):
            break;
        default:
            GMX_THROW(InternalError("Enum for frameTime not recognized"));
    }

    haveProcessedFirstFrame_ = true;
}

void
SetTime::setTimeShift(real newFrameTime)
{
    // get difference between frames for new time.
    if (frameTime_ == ChangeFrameTimeType::efStartTime)
    {
        timeShift_ = newFrameTime - timeStep_;
    }
}

} // namespace gmx
