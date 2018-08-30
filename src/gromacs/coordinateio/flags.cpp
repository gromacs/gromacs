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
 * Implements gmx::CoordinateFileWriteFlags.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "flags.h"

#include <algorithm>

#include "gromacs/compat/make_unique.h"
#include "gromacs/coordinateio/modules/setatoms.h"
#include "gromacs/coordinateio/modules/setbox.h"
#include "gromacs/coordinateio/modules/setforces.h"
#include "gromacs/coordinateio/modules/setprecision.h"
#include "gromacs/coordinateio/modules/settime.h"
#include "gromacs/coordinateio/modules/setvelocities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"

namespace gmx
{

void
CoordinateFileWriteFlags::initFileOptions(IOptionsContainer *options)
{
    options->addOption(EnumOption<ChangeSettingType>("vel")
                           .enumValue(cChangeSettingTypeEnum)
                           .store(&velocity_)
                           .description("Save velocities from frame if possible"));
    options->addOption(EnumOption<ChangeSettingType>("force")
                           .enumValue(cChangeSettingTypeEnum)
                           .store(&force_)
                           .description("Save forces from frame if possible"));
    options->addOption(EnumOption<ChangeSettingType>("atoms")
                           .enumValue(cChangeSettingTypeEnum)
                           .store(&atoms_)
                           .description("Decide on providing new atom information from topology or using current frame atom information"));
    options->addOption(EnumOption<ChangeFrameUnchangedYesType>("precision")
                           .enumValue(cChangeFrameUnchangedYesTypeEnum)
                           .store(&precision_)
                           .description("Decide on changing current frame precision or not"));
    options->addOption(IntegerOption("newprec")
                           .store(&prec_)
                           .defaultValue(3)
                           .description("Set output precision to custom value"));
    options->addOption(EnumOption<ChangeFrameTimeType>("frametime")
                           .enumValue(cChangeFrameTimeTypeEnum)
                           .store(&frameTime_)
                           .description("Decide on if time information in a coordinate frame should changed or not"));
    options->addOption(DoubleOption("starttime")
                           .store(&startTimeValue_)
                           .defaultValue(0)
                           .timeValue()
                           .description("Change start time for first frame"));
    options->addOption(DoubleOption("timestep")
                           .store(&timeStepValue_)
                           .defaultValue(0)
                           .timeValue()
                           .description("Change time between different frames"));
    options->addOption(EnumOption<ChangeFrameUnchangedYesType>("box")
                           .enumValue(cChangeFrameUnchangedYesTypeEnum)
                           .store(&box_)
                           .description("Decide on if frame box information should be changed"));
    options->addOption(DoubleOption("newbox")
                           .vector()
                           .storeVector(&newBoxVector_)
                           .valueCount(3)
                           .description("New diagonal box vector for output frame"));
}

void
CoordinateFileWriteFlags::checkOptions()
{
    if (box_ != ChangeFrameUnchangedYesType::efUnchanged)
    {
        for (int i = 0; i < DIM; i++)
        {
            newBox_[i][i] = newBoxVector_[i];
        }
    }
}

void
CoordinateFileWriteFlags::registerModules(const OutputManagerPointer &output)
{
    if (velocity_ != ChangeSettingType::efUnchanged)
    {
        output->addFlagModule(compat::make_unique<SetVelocities>(velocity_));
    }
    if (force_ != ChangeSettingType::efUnchanged)
    {
        output->addFlagModule(compat::make_unique<SetForces>(force_));
    }
    if (precision_ != ChangeFrameUnchangedYesType::efUnchanged)
    {
        output->addFlagModule(compat::make_unique<SetPrecision>(prec_));
    }
    if (atoms_ != ChangeSettingType::efUnchanged)
    {
        output->addFlagModule(compat::make_unique<SetAtoms>(atoms_));
    }
    if (frameTime_ != ChangeFrameTimeType::efUnchanged)
    {
        output->addFlagModule(compat::make_unique<SetTime>(startTimeValue_, timeStepValue_, frameTime_));
    }
    if (box_ != ChangeFrameUnchangedYesType::efUnchanged)
    {
        output->addFlagModule(compat::make_unique<SetBox>(newBox_));
    }
}


} // namespace gmx
