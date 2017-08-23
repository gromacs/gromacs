/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests for NonbondedOnGpuFromUser
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/nonbondedongpufromuser.h"

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/refdata.h"

#include "nodetaskassigner-common.h"

namespace gmx
{

namespace test
{
namespace
{

using testing::SizeIs;

TEST(NonbondedOnGpuFromUserChooserTest, ThrowsWhenGpusWontWork)
{
    for (const auto &gpuIdString : allUserAssignedGpuDeviceIds)
    {
        for (const auto &gpuCompatibilitiesPair : allGpuCompatibilities)
        {
            gmx_hw_info_t  hardwareInfo;
            hardwareInfo.compatibleGpus = gpuCompatibilitiesPair.second;
            auto           taskAssigner   = compat::make_unique<NonbondedOnGpuFromUser>(gpuIdString, hardwareInfo);

            // To use designated GPUs for nonbondeds, the user must
            // choose a calculation type where that is supported.
            EXPECT_THROW(taskAssigner->decideWhetherToUseGpus(true,  false), InconsistentInputError);
            EXPECT_THROW(taskAssigner->decideWhetherToUseGpus(false, true),  InconsistentInputError);
            EXPECT_THROW(taskAssigner->decideWhetherToUseGpus(false, false), InconsistentInputError);
        }
    }
}

//! Parameterized test fixture for when an automated task assignment is made
class NonbondedOnGpuFromUserTest : public ::testing::Test
{
    public:
        //! Test reference data.
        TestReferenceData    refData_;
        //! Test reference data checker.
        TestReferenceChecker rootChecker_;
        //! Constructor.
        NonbondedOnGpuFromUserTest() : refData_(),
                                       rootChecker_(refData_.rootChecker()) {}
        //! Helper for running the tests.
        void runTest(NonbondedOnGpuFromUser *taskAssigner,
                     TestReferenceChecker   *checker);
};

void NonbondedOnGpuFromUserTest::runTest(NonbondedOnGpuFromUser *taskAssigner,
                                         TestReferenceChecker   *checker)
{
    EXPECT_NO_THROW(taskAssigner->decideWhetherToUseGpus(true, true));
    EXPECT_TRUE(taskAssigner->areNonbondedOnGpu());

    for (const auto &taskEligibility : allTaskEligibilities)
    {
        auto               gpuTasksOnRanks = taskEligibility.second;
        GpuTaskAssignments assignments;
        try
        {
            assignments = taskAssigner->assignGpuTasksToDeviceIds(gpuTasksOnRanks);
        }
        catch (InconsistentInputError &ex)
        {
            // There is no further behaviour that we can test, given invalid inputs.
            continue;
        }

        auto compound(checker->checkCompound("AssignmentResults", nullptr));
        compound.checkString(taskEligibility.first, "Given Gpu Tasks On Ranks");
        EXPECT_THAT(assignments, SizeIs(gpuTasksOnRanks.size())) << "Valid assignment on each rank";
        compound.checkSequence(assignments.begin(), assignments.end(), "Produced Assignments On Ranks", gpuTaskAssignmentsItemChecker);

        EXPECT_TRUE(taskAssigner->didUserAssignGpus());
        // TODO There could be log output that we might assert upon
    }
}

TEST_F(NonbondedOnGpuFromUserTest, MakesAssignmentsAsExpected)
{
    // Note that the reference data for these tests are sensitive to
    // the order in the containers of test inputs, so that when
    // extending or changing them, then the reference data might need
    // to be regenerated.
    for (const auto &gpuIdString : allUserAssignedGpuDeviceIds)
    {
        for (const auto &gpuCompatibilitiesPair : allGpuCompatibilities)
        {
            gmx_hw_info_t        hardwareInfo;
            hardwareInfo.compatibleGpus = gpuCompatibilitiesPair.second;
            auto                 taskAssigner   = compat::make_unique<NonbondedOnGpuFromUser>(gpuIdString, hardwareInfo);

            TestReferenceChecker checker(rootChecker_.checkCompound("InputCase", nullptr));
            checker.checkString(gpuIdString, "GpuIdString");
            checker.checkString(gpuCompatibilitiesPair.first, "Given GPU compatibilities");
            runTest(taskAssigner.get(), &checker);
        }
    }
}

} // namespace
} // namespace
} // namespace
