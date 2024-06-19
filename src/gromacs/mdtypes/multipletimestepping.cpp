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
#include "gmxpre.h"

#include "multipletimestepping.h"

#include <memory>
#include <optional>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

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

std::vector<MtsLevel> setupMtsLevels(const GromppMtsOpts& mtsOpts, std::vector<std::string>* errorMessages)
{
    std::vector<MtsLevel> mtsLevels;

    if (mtsOpts.numLevels != 2)
    {
        if (errorMessages)
        {
            errorMessages->emplace_back("Only mts-levels = 2 is supported");
        }
    }
    else
    {
        mtsLevels.resize(2);

        const std::vector<std::string> inputForceGroups = gmx::splitString(mtsOpts.level2Forces);
        auto&                          forceGroups      = mtsLevels[1].forceGroups;
        for (const auto& inputForceGroup : inputForceGroups)
        {
            bool found     = false;
            int  nameIndex = 0;
            for (const auto& forceGroupName : gmx::mtsForceGroupNames)
            {
                if (gmx::equalCaseInsensitive(inputForceGroup, forceGroupName))
                {
                    forceGroups.set(nameIndex);
                    found = true;
                }
                nameIndex++;
            }
            if (!found && errorMessages)
            {
                errorMessages->push_back(
                        gmx::formatString("Unknown MTS force group '%s'", inputForceGroup.c_str()));
            }
        }

        // Make the level 0 use the complement of the force groups of group 1
        mtsLevels[0].forceGroups = ~mtsLevels[1].forceGroups;
        mtsLevels[0].stepFactor  = 1;

        mtsLevels[1].stepFactor = mtsOpts.level2Factor;

        if (errorMessages && mtsLevels[1].stepFactor <= 1)
        {
            errorMessages->emplace_back("mts-factor should be larger than 1");
        }
    }

    return mtsLevels;
}

bool haveValidMtsSetup(const t_inputrec& ir)
{
    return (ir.useMts && ir.mtsLevels.size() == 2 && ir.mtsLevels[1].stepFactor > 1);
}

namespace
{

//! Checks whether \p nstValue is a multiple of the largest MTS step, returns an error string for parameter \p param when this is not the case
std::optional<std::string> checkMtsInterval(ArrayRef<const MtsLevel> mtsLevels, const char* param, const int nstValue)
{
    GMX_RELEASE_ASSERT(mtsLevels.size() >= 2, "Need at least two levels for MTS");

    const int mtsFactor = mtsLevels.back().stepFactor;
    if (nstValue % mtsFactor == 0)
    {
        return {};
    }
    else
    {
        return gmx::formatString(
                "With MTS, %s = %d should be a multiple of mts-factor = %d", param, nstValue, mtsFactor);
    }
}

} // namespace

std::vector<std::string> checkMtsRequirements(const t_inputrec& ir)
{
    std::vector<std::string> errorMessages;

    if (!ir.useMts)
    {
        return errorMessages;
    }

    GMX_RELEASE_ASSERT(haveValidMtsSetup(ir), "MTS setup should be valid here");

    ArrayRef<const MtsLevel> mtsLevels = ir.mtsLevels;

    if (ir.eI != IntegrationAlgorithm::MD)
    {
        errorMessages.push_back(
                gmx::formatString("Multiple time stepping is only supported with integrator %s",
                                  enumValueToString(IntegrationAlgorithm::MD)));
    }

    if ((usingFullElectrostatics(ir.coulombtype) || usingLJPme(ir.vdwtype))
        && forceGroupMtsLevel(ir.mtsLevels, MtsForceGroups::LongrangeNonbonded) == 0)
    {
        errorMessages.emplace_back(
                "With long-range electrostatics and/or LJ treatment, the long-range part "
                "has to be part of the mts-level2-forces");
    }

    std::optional<std::string> mesg;
    if (ir.nstcalcenergy > 0)
    {
        if ((mesg = checkMtsInterval(mtsLevels, "nstcalcenergy", ir.nstcalcenergy)))
        {
            errorMessages.push_back(mesg.value());
        }
    }
    if ((mesg = checkMtsInterval(mtsLevels, "nstenergy", ir.nstenergy)))
    {
        errorMessages.push_back(mesg.value());
    }
    if ((mesg = checkMtsInterval(mtsLevels, "nstlog", ir.nstlog)))
    {
        errorMessages.push_back(mesg.value());
    }
    if ((mesg = checkMtsInterval(mtsLevels, "nstfout", ir.nstfout)))
    {
        errorMessages.push_back(mesg.value());
    }
    if (ir.efep != FreeEnergyPerturbationType::No)
    {
        if ((mesg = checkMtsInterval(mtsLevels, "nstdhdl", ir.fepvals->nstdhdl)))
        {
            errorMessages.push_back(mesg.value());
        }
    }
    if (mtsLevels.back().forceGroups[static_cast<int>(gmx::MtsForceGroups::Nonbonded)])
    {
        if ((mesg = checkMtsInterval(mtsLevels, "nstlist", ir.nstlist)))
        {
            errorMessages.push_back(mesg.value());
        }
    }

    if (ir.bPull)
    {
        const int pullMtsLevel  = gmx::forceGroupMtsLevel(ir.mtsLevels, gmx::MtsForceGroups::Pull);
        const int mtsStepFactor = ir.mtsLevels[pullMtsLevel].stepFactor;
        if (ir.pull->nstxout % mtsStepFactor != 0)
        {
            errorMessages.emplace_back("pull-nstxout should be a multiple of mts-factor");
        }
        if (ir.pull->nstfout % mtsStepFactor != 0)
        {
            errorMessages.emplace_back("pull-nstfout should be a multiple of mts-factor");
        }
    }

    return errorMessages;
}

} // namespace gmx
