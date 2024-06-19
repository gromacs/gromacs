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
/*!\internal
 * \file
 * \brief
 * Implements helper function to populate requirements from user input.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "requirements.h"

#include <algorithm>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

//! Mapping for enums from \ref ChangeSettingType.
static const EnumerationArray<ChangeSettingType, const char*> c_changeSettingTypeNames = {
    { "preserved-if-present", "always", "never" }
};
//! Mapping for enums from \ref ChangeAtomsType.
static const EnumerationArray<ChangeAtomsType, const char*> c_changeAtomsTypeNames = {
    { "preserved-if-present", "always-from-structure", "never", "always" }
};
/* Currently unused
//! Mapping for enums from \ref ChangeFrameInfoType.
static const EnumerationArray<ChangeFrameInfoType, const char*> c_changeFrameInfoTypeNames = { {
"preserved-if-present", "always" } };
//! Mapping for values from \ref ChangeFrameTimeType.
static const EnumerationArray<ChangeFrameTimeType, const char*> c_changeFrameTimeTypeNames = { {
"preserved-if-present", "starttime", "timestep", "both" } };
*/

void OutputRequirementOptionDirector::initOptions(IOptionsContainer* options)
{
    options->addOption(EnumOption<ChangeSettingType>("vel")
                               .enumValue(c_changeSettingTypeNames)
                               .store(&velocity_)
                               .description("Save velocities from frame if possible"));
    options->addOption(EnumOption<ChangeSettingType>("force")
                               .enumValue(c_changeSettingTypeNames)
                               .store(&force_)
                               .description("Save forces from frame if possible"));
    options->addOption(
            EnumOption<ChangeAtomsType>("atoms").enumValue(c_changeAtomsTypeNames).store(&atoms_).description("Decide on providing new atom information from topology or using current frame atom information"));
    options->addOption(IntegerOption("precision")
                               .store(&prec_)
                               .defaultValue(prec_)
                               .storeIsSet(&setNewPrecision_)
                               .description("Set output precision to custom value"));
    options->addOption(RealOption("starttime")
                               .store(&startTimeValue_)
                               .defaultValue(startTimeValue_)
                               .timeValue()
                               .storeIsSet(&setNewStartTime_)
                               .description("Change start time for first frame"));
    options->addOption(RealOption("timestep")
                               .store(&timeStepValue_)
                               .defaultValue(timeStepValue_)
                               .timeValue()
                               .storeIsSet(&setNewTimeStep_)
                               .description("Change time between different frames"));
    options->addOption(RealOption("box")
                               .vector()
                               .storeVector(&newBoxVector_)
                               .valueCount(3)
                               .storeIsSet(&setNewBox_)
                               .description("New diagonal box vector for output frame"));
}

OutputRequirements OutputRequirementOptionDirector::process() const
{
    OutputRequirements requirements;
    /* If the user has just set the values directly without setting the flags,
     * we set the flags to state that user requested changes are there.*/
    if (setNewBox_)
    {
        requirements.box = ChangeFrameInfoType::Always;
        clear_mat(requirements.newBox);
        for (int i = 0; i < DIM; ++i)
        {
            requirements.newBox[i][i] = newBoxVector_[i];
        }
    }
    if (setNewPrecision_)
    {
        requirements.precision = ChangeFrameInfoType::Always;
        requirements.prec      = prec_;
    }
    if ((setNewTimeStep_ || setNewStartTime_))
    {
        requirements.startTimeValue = startTimeValue_;
        requirements.timeStepValue  = timeStepValue_;
        if (setNewTimeStep_ && setNewStartTime_)
        {
            requirements.frameTime = ChangeFrameTimeType::Both;
        }
        else if (setNewTimeStep_)
        {
            requirements.frameTime = ChangeFrameTimeType::TimeStep;
        }
        else
        {
            requirements.frameTime = ChangeFrameTimeType::StartTime;
        }
    }
    requirements.atoms    = atoms_;
    requirements.velocity = velocity_;
    requirements.force    = force_;
    return requirements;
}

} // namespace gmx
