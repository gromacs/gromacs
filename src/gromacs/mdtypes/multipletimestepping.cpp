/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "multipletimestepping.h"

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

int nonbondedMtsFactor(const t_inputrec& ir)
{
    GMX_RELEASE_ASSERT(!ir.useMts || ir.mtsLevels.size() == 2, "Only 2 MTS levels supported here");

    if (ir.useMts && ir.mtsLevels[1].forceGroups[static_cast<int>(MtsForceGroups::Nonbonded)])
    {
        return ir.mtsLevels[1].stepFactor;
    }
    else
    {
        return 1;
    }
}

void assertMtsRequirements(const t_inputrec& ir)
{
    if (!ir.useMts)
    {
        return;
    }

    GMX_RELEASE_ASSERT(ir.mtsLevels.size() >= 2, "Need at least two levels for MTS");

    GMX_RELEASE_ASSERT(ir.mtsLevels[0].stepFactor == 1, "Base MTS step should be 1");

    GMX_RELEASE_ASSERT(
            (!EEL_FULL(ir.coulombtype) && !EVDW_PME(ir.vdwtype))
                    || ir.mtsLevels.back().forceGroups[static_cast<int>(MtsForceGroups::LongrangeNonbonded)],
            "Long-range nonbondeds should be in the highest MTS level");

    for (const auto& mtsLevel : ir.mtsLevels)
    {
        const int mtsFactor = mtsLevel.stepFactor;
        GMX_RELEASE_ASSERT(ir.nstcalcenergy % mtsFactor == 0,
                           "nstcalcenergy should be a multiple of mtsFactor");
        GMX_RELEASE_ASSERT(ir.nstenergy % mtsFactor == 0,
                           "nstenergy should be a multiple of mtsFactor");
        GMX_RELEASE_ASSERT(ir.nstlog % mtsFactor == 0, "nstlog should be a multiple of mtsFactor");
        GMX_RELEASE_ASSERT(ir.epc == epcNO || ir.nstpcouple % mtsFactor == 0,
                           "nstpcouple should be a multiple of mtsFactor");
        GMX_RELEASE_ASSERT(ir.efep == efepNO || ir.fepvals->nstdhdl % mtsFactor == 0,
                           "nstdhdl should be a multiple of mtsFactor");
        if (ir.mtsLevels.back().forceGroups[static_cast<int>(gmx::MtsForceGroups::Nonbonded)])
        {
            GMX_RELEASE_ASSERT(ir.nstlist % ir.mtsLevels.back().stepFactor == 0,
                               "With multiple time stepping for the non-bonded pair interactions, "
                               "nstlist should be a "
                               "multiple of mtsFactor");
        }
    }
}

} // namespace gmx
