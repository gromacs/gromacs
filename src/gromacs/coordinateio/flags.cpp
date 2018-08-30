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
 * Implements helper function to populate requirements from user input.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "flags.h"

#include <algorithm>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

void
initOutputRequirementFileOptions(IOptionsContainer  *options,
                                 OutputRequirements *requirements)
{
    options->addOption(EnumOption<ChangeSettingType>("vel")
                           .enumValue(cChangeSettingTypeEnum)
                           .store(&requirements->velocity)
                           .description("Save velocities from frame if possible"));
    options->addOption(EnumOption<ChangeSettingType>("force")
                           .enumValue(cChangeSettingTypeEnum)
                           .store(&requirements->force)
                           .description("Save forces from frame if possible"));
    options->addOption(EnumOption<ChangeAtomsType>("atoms")
                           .enumValue(cChangeAtomsTypeEnum)
                           .store(&requirements->atoms)
                           .description("Decide on providing new atom information from topology or using current frame atom information"));
    options->addOption(EnumOption<ChangeFrameInfoType>("precision")
                           .enumValue(cChangeFrameInfoTypeEnum)
                           .store(&requirements->precision)
                           .description("Decide on changing current frame precision or not"));
    options->addOption(IntegerOption("newprec")
                           .store(&requirements->prec)
                           .defaultValue(requirements->prec)
                           .storeIsSet(&requirements->setNewPrecision)
                           .description("Set output precision to custom value"));
    options->addOption(EnumOption<ChangeFrameTimeType>("frametime")
                           .enumValue(cChangeFrameTimeTypeEnum)
                           .store(&requirements->frameTime)
                           .description("Decide on if time information in a coordinate frame should changed or not"));
    options->addOption(RealOption("starttime")
                           .store(&requirements->startTimeValue)
                           .defaultValue(requirements->startTimeValue)
                           .timeValue()
                           .storeIsSet(&requirements->setNewStartTime)
                           .description("Change start time for first frame"));
    options->addOption(RealOption("timestep")
                           .store(&requirements->timeStepValue)
                           .defaultValue(requirements->timeStepValue)
                           .timeValue()
                           .storeIsSet(&requirements->setNewTimeStep)
                           .description("Change time between different frames"));
    options->addOption(EnumOption<ChangeFrameInfoType>("box")
                           .enumValue(cChangeFrameInfoTypeEnum)
                           .store(&requirements->box)
                           .description("Decide on if frame box information should be changed"));
    options->addOption(RealOption("newbox")
                           .vector()
                           .storeVector(&requirements->newBoxVector)
                           .valueCount(3)
                           .storeIsSet(&requirements->setNewBox)
                           .description("New diagonal box vector for output frame"));
}

void
checkRequirementsOptions(OutputRequirements *reqs)
{
    /* If the user has just set the values directly without setting the flags,
     * we set the flags to state that user requested changes are there.*/
    if (reqs->setNewBox && reqs->box == ChangeFrameInfoType::efUnchanged)
    {
        reqs->box = ChangeFrameInfoType::efUserYes;
    }
    if (reqs->setNewPrecision && reqs->precision == ChangeFrameInfoType::efUnchanged)
    {
        reqs->precision = ChangeFrameInfoType::efUserYes;
    }
    if ((reqs->setNewTimeStep || reqs->setNewStartTime) &&
        reqs->frameTime == ChangeFrameTimeType::efUnchanged)
    {
        if (reqs->setNewTimeStep && reqs->setNewStartTime)
        {
            reqs->frameTime = ChangeFrameTimeType::efBothTime;
        }
        else if (reqs->setNewTimeStep)
        {
            reqs->frameTime = ChangeFrameTimeType::efTimeStep;
        }
        else
        {
            reqs->frameTime = ChangeFrameTimeType::efStartTime;
        }
    }
    if (reqs->box != ChangeFrameInfoType::efUnchanged)
    {
        if (!reqs->setNewBox)
        {
            GMX_THROW(InconsistentInputError("Need box vector to assign new box"));
        }
        for (int i = 0; i < DIM; i++)
        {
            reqs->newBox[i][i] = reqs->newBoxVector[i];
        }
    }
    reqs->isValid = true;
}

} // namespace gmx
