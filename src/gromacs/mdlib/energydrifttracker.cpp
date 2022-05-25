/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Implements functions from the EnergyDriftTracker class.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "energydrifttracker.h"

#include <cmath>

#include <string>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

void EnergyDriftTracker::addPoint(double time, double energy)
{
    GMX_ASSERT(std::isfinite(energy), "Non-finite energy encountered!");

    if (!storedFirst_)
    {
        firstTime_   = time;
        firstEnergy_ = energy;
        storedFirst_ = true;
    }
    lastTime_   = time;
    lastEnergy_ = energy;
}

double EnergyDriftTracker::energyDrift() const
{
    if (timeInterval() > 0)
    {
        return (lastEnergy_ - firstEnergy_) / (timeInterval() * numAtoms_);
    }
    else
    {
        return 0;
    }
}

std::string EnergyDriftTracker::energyDriftString(const std::string& partName) const
{
    std::string mesg;

    if (timeInterval() > 0)
    {
        mesg = formatString("Energy conservation over %s of length %g ps, time %g to %g ps\n",
                            partName.c_str(),
                            timeInterval(),
                            firstTime_,
                            lastTime_);
        mesg += formatString("  Conserved energy drift: %.2e kJ/mol/ps per atom\n", energyDrift());
    }
    else
    {
        mesg = formatString(
                "Time interval for measuring conserved energy has length 0, time %g to %g ps\n",
                firstTime_,
                lastTime_);
    }

    return mesg;
}

} // namespace gmx
