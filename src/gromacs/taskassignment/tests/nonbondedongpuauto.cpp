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
 * Tests for NonbondedOnGpuAuto.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/nonbondedongpuauto.h"

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

using testing::IsEmpty;
using testing::SizeIs;

//! Helper array for range-based for loops.
const bool booleans_g[] = { false, true };

TEST(NonbondedOnGpuAutoChooserTest, ThrowsWhenGpusWontWork)
{
    for (const auto &gpuCompatibilitiesPair : allGpuCompatibilities)
    {
        gmx_hw_info_t hardwareInfo;
        hardwareInfo.compatibleGpus = gpuCompatibilitiesPair.second;
        {
            bool forceNonbondedOnPhysicalGpu = true;
            auto taskAssigner                = compat::make_unique<NonbondedOnGpuAuto>(forceNonbondedOnPhysicalGpu, hardwareInfo);

            // When forcing nonbondeds to run on GPUs, the user
            // must choose a calculation type where that is
            // supported.
            EXPECT_THROW(taskAssigner->chooseWhetherToUseGpus(true,  false), InconsistentInputError);
            EXPECT_THROW(taskAssigner->chooseWhetherToUseGpus(false, true),  InconsistentInputError);
            EXPECT_THROW(taskAssigner->chooseWhetherToUseGpus(false, false), InconsistentInputError);
        }
        {
            bool forceNonbondedOnPhysicalGpu = false;
            auto taskAssigner                = compat::make_unique<NonbondedOnGpuAuto>(forceNonbondedOnPhysicalGpu, hardwareInfo);

            EXPECT_NO_THROW(taskAssigner->chooseWhetherToUseGpus(true,  false));
            EXPECT_FALSE(taskAssigner->areNonbondedOnGpu());

            EXPECT_NO_THROW(taskAssigner->chooseWhetherToUseGpus(false, true));
            EXPECT_FALSE(taskAssigner->areNonbondedOnGpu());

            EXPECT_NO_THROW(taskAssigner->chooseWhetherToUseGpus(false, false));
            EXPECT_FALSE(taskAssigner->areNonbondedOnGpu());
        }
    }
}

//! Parameterized test fixture for when an automated task assignment is made
class NonbondedOnGpuAutoTest : public ::testing::Test
{
    public:
        //! Test reference data.
        TestReferenceData    refData_;
        //! Test reference data checker.
        TestReferenceChecker rootChecker_;
        //! Constructor.
        NonbondedOnGpuAutoTest() : refData_(),
                                   rootChecker_(refData_.rootChecker()) {}
        //! Helper for running the tests.
        void runTest(bool forceNonbondedOnPhysicalGpu);
};

void NonbondedOnGpuAutoTest::runTest(bool forceNonbondedOnPhysicalGpu)
{
    // Note that the reference data for these tests are sensitive to
    // the order in the containers of test inputs, so that when
    // extending or changing them, then the reference data might need
    // to be regenerated.
    for (const auto &gpuCompatibilitiesPair : allGpuCompatibilities)
    {
        gmx_hw_info_t hardwareInfo;
        hardwareInfo.compatibleGpus = gpuCompatibilitiesPair.second;
        auto          taskAssigner   = compat::make_unique<NonbondedOnGpuAuto>(forceNonbondedOnPhysicalGpu, hardwareInfo);

        for (const auto &usingVerletScheme : booleans_g)
        {
            for (const auto &gpuAccelerationIsUseful : booleans_g)
            {
                try
                {
                    taskAssigner->chooseWhetherToUseGpus(usingVerletScheme, gpuAccelerationIsUseful);
                }
                catch (InconsistentInputError &ex)
                {
                    if (forceNonbondedOnPhysicalGpu)
                    {
                        // This is an expected behaviour. When forcing
                        // nonbondeds to run on GPUs, the user must
                        // choose a calculation type where that is
                        // supported. There is no further behaviour
                        // that we can test, given these inconsistent
                        // inputs.
                        continue;
                    }
                    else
                    {
                        ADD_FAILURE() << "Unexpected throw from chooseWhetherToUseGpus";
                    }
                }

                // Note that assignGpuTasksToDeviceIds is const, so
                // it's ok to keep re-using the task assigner for
                // different sets of eligible tasks.
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
                        // The available tasks don't match the devices, so
                        // there is no further behaviour that we can test,
                        // given these inconsistent inputs.
                        continue;
                    }

                    TestReferenceChecker checker(rootChecker_.checkCompound("InputCase", nullptr));
                    checker.checkBoolean(forceNonbondedOnPhysicalGpu, "ForceNonbondedOnPhysicalGpu");
                    checker.checkString(gpuCompatibilitiesPair.first, "Given GPU compatibilities");
                    checker.checkString(taskEligibility.first, "Given Gpu Tasks On Ranks");

                    EXPECT_THAT(assignments, SizeIs(gpuTasksOnRanks.size())) << "Valid assignment on each rank";
                    checker.checkSequence(assignments.begin(), assignments.end(), "Produced Assignments On Ranks", gpuTaskAssignmentsItemChecker);

                    EXPECT_FALSE(taskAssigner->didUserAssignGpus());
                    // TODO There could be log output that we might assert upon
                }
            }
        }
    }
}

TEST_F(NonbondedOnGpuAutoTest, MakesAssignmentsAsExpected)
{
    {
        bool forceNonbondedOnPhysicalGpu = true;
        runTest(forceNonbondedOnPhysicalGpu);
    }
    {
        bool forceNonbondedOnPhysicalGpu = false;
        runTest(forceNonbondedOnPhysicalGpu);
    }
}

} // namespace
} // namespace
} // namespace
