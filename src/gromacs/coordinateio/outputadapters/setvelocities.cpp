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
 * Implements SetVelocities class.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "setvelocities.h"

#include <string>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

void SetVelocities::checkAbilityDependencies(unsigned long abilities) const
{
    if ((abilities & convertFlag(moduleRequirements_)) == 0U)
    {
        std::string errorMessage =
                "Output file type does not support writing velocities. "
                "Only GRO, TRR and TNG support this output.";
        GMX_THROW(InconsistentInputError(errorMessage.c_str()));
    }
}

void SetVelocities::processFrame(const int /*framenumber*/, t_trxframe* input)
{
    switch (velocity_)
    {
        case (ChangeSettingType::Never):
            input->bV = false;
            input->v  = nullptr;
            break;
        case (ChangeSettingType::Always):
            if (!input->bV)
            {
                GMX_THROW(InconsistentInputError(
                        "Velocity output requested but current frame has no velocities"));
            }
            break;
        default: GMX_THROW(InconsistentInputError("Value for velocity flag is not supported"));
    }
}

} // namespace gmx
