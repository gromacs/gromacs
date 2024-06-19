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
 * \brief
 * Tests for the MultipleTimeStepping class and stand-alone functions.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/multipletimestepping.h"

#include <bitset>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

//! brief Sets up the MTS levels in \p ir and tests whether the number of errors matches \p numExpectedErrors
void setAndCheckMtsLevels(const GromppMtsOpts& mtsOpts, t_inputrec* ir, const int numExpectedErrors)
{
    std::vector<std::string> errorMessages;
    ir->useMts    = true;
    ir->mtsLevels = setupMtsLevels(mtsOpts, &errorMessages);

    if (haveValidMtsSetup(*ir))
    {
        std::vector<std::string> errorMessagesCheck = checkMtsRequirements(*ir);

        // Concatenate the two lists with error messages
        errorMessages.insert(errorMessages.end(), errorMessagesCheck.begin(), errorMessagesCheck.end());
    }

    EXPECT_EQ(errorMessages.size(), numExpectedErrors);
}

} // namespace

//! Checks that only numLevels = 2 does not produce an error
TEST(MultipleTimeStepping, ChecksNumLevels)
{
    for (int numLevels = 0; numLevels <= 3; numLevels++)
    {
        GromppMtsOpts mtsOpts;
        mtsOpts.numLevels    = numLevels;
        mtsOpts.level2Factor = 2;

        t_inputrec ir;

        setAndCheckMtsLevels(mtsOpts, &ir, numLevels != 2 ? 1 : 0);
    }
}

//! Test that each force group works
TEST(MultipleTimeStepping, SelectsForceGroups)
{
    for (int forceGroupIndex = 0; forceGroupIndex < static_cast<int>(MtsForceGroups::Count);
         forceGroupIndex++)
    {
        const MtsForceGroups forceGroup = static_cast<MtsForceGroups>(forceGroupIndex);
        SCOPED_TRACE("Testing force group " + mtsForceGroupNames[forceGroup]);

        GromppMtsOpts mtsOpts;
        mtsOpts.numLevels    = 2;
        mtsOpts.level2Forces = mtsForceGroupNames[forceGroup];
        mtsOpts.level2Factor = 2;

        t_inputrec ir;

        setAndCheckMtsLevels(mtsOpts, &ir, 0);

        EXPECT_EQ(ir.mtsLevels[1].forceGroups.count(), 1);
        EXPECT_EQ(ir.mtsLevels[1].forceGroups[forceGroupIndex], true);
    }
}

//! Checks that factor is checked
TEST(MultipleTimeStepping, ChecksStepFactor)
{
    for (int stepFactor = 0; stepFactor <= 3; stepFactor++)
    {
        GromppMtsOpts mtsOpts;
        mtsOpts.numLevels    = 2;
        mtsOpts.level2Factor = stepFactor;

        t_inputrec ir;

        setAndCheckMtsLevels(mtsOpts, &ir, stepFactor < 2 ? 1 : 0);
    }
}

namespace
{

GromppMtsOpts simpleMtsOpts()
{
    GromppMtsOpts mtsOpts;
    mtsOpts.numLevels    = 2;
    mtsOpts.level2Forces = "nonbonded";
    mtsOpts.level2Factor = 4;

    return mtsOpts;
}

} // namespace

TEST(MultipleTimeStepping, ChecksPmeIsAtLastLevel)
{
    const GromppMtsOpts mtsOpts = simpleMtsOpts();

    t_inputrec ir;
    ir.coulombtype = CoulombInteractionType::Pme;

    setAndCheckMtsLevels(mtsOpts, &ir, 1);
}

//! Test fixture base for parametrizing interval tests
using MtsIntervalTestParams = std::tuple<std::string, int>;
class MtsIntervalTest : public ::testing::Test, public ::testing::WithParamInterface<MtsIntervalTestParams>
{
public:
    MtsIntervalTest()
    {
        const auto  params        = GetParam();
        const auto& parameterName = std::get<0>(params);
        const auto  interval      = std::get<1>(params);
        numExpectedErrors_        = (interval == 4 ? 0 : 1);

        if (parameterName == "nstcalcenergy")
        {
            ir_.nstcalcenergy = interval;
        }
        else if (parameterName == "nstenergy")
        {
            ir_.nstenergy = interval;
        }
        else if (parameterName == "nstfout")
        {
            ir_.nstfout = interval;
        }
        else if (parameterName == "nstlist")
        {
            ir_.nstlist = interval;
        }
        else if (parameterName == "nstdhdl")
        {
            ir_.efep             = FreeEnergyPerturbationType::Yes;
            ir_.fepvals->nstdhdl = interval;
        }
        else

        {
            GMX_RELEASE_ASSERT(false, "unknown parameter name");
        }
    }

    t_inputrec ir_;
    int        numExpectedErrors_;
};

TEST_P(MtsIntervalTest, Works)
{
    const GromppMtsOpts mtsOpts = simpleMtsOpts();

    setAndCheckMtsLevels(mtsOpts, &ir_, numExpectedErrors_);
}

INSTANTIATE_TEST_SUITE_P(
        ChecksStepInterval,
        MtsIntervalTest,
        ::testing::Combine(
                ::testing::Values("nstcalcenergy", "nstenergy", "nstfout", "nstlist", "nstdhdl"),
                ::testing::Values(3, 4, 5)));

// Check that correct input does not produce errors
TEST(MultipleTimeStepping, ChecksIntegrator)
{
    const GromppMtsOpts mtsOpts = simpleMtsOpts();

    t_inputrec ir;
    ir.eI = IntegrationAlgorithm::BD;

    setAndCheckMtsLevels(mtsOpts, &ir, 1);
}

} // namespace test
} // namespace gmx
