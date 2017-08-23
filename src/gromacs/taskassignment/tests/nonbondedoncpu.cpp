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
 * Tests for NonbondedOnCpu.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/nonbondedoncpu.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"

#include "nodetaskassigner-common.h"

namespace gmx
{

namespace test
{
namespace
{

using testing::Bool;
using testing::Combine;
using testing::IsEmpty;
using testing::SizeIs;
using testing::Values;
using testing::ValuesIn;

//! Convenience typedef for parameterized testing.
using TestParameters = std::tuple<bool,
                                  bool>;

//! Parameterized test fixture for when a task assignment is made that uses no GPUs.
class NonbondedOnCpuTest : public ::testing::TestWithParam<TestParameters>
{
    public:
        //! Argument passed to the assigner.
        bool usingVerletScheme_;
        //! Argument passed to the assigner.
        bool gpuAccelerationIsUseful_;
        //! The task assigner.
        std::unique_ptr<NonbondedOnCpu> taskAssigner_;
        //! Constructor.
        NonbondedOnCpuTest() : usingVerletScheme_(std::get<0>(GetParam())),
                               gpuAccelerationIsUseful_(std::get<1>(GetParam())),
                               taskAssigner_(compat::make_unique<NonbondedOnCpu>()) {}
};

TEST_P(NonbondedOnCpuTest, MakesEmptyAssignment)
{
    EXPECT_NO_THROW(taskAssigner_->chooseWhetherToUseGpus(usingVerletScheme_, gpuAccelerationIsUseful_));
    EXPECT_FALSE(taskAssigner_->areNonbondedOnGpu());

    for (const auto &taskEligibility : allTaskEligibilities)
    {
        auto               gpuTasksOnRanks = taskEligibility.second;
        GpuTaskAssignments assignments;
        EXPECT_NO_THROW(assignments = taskAssigner_->assignGpuTasksToDeviceIds(gpuTasksOnRanks));
        EXPECT_THAT(assignments, SizeIs(gpuTasksOnRanks.size())) << "Valid assignment on each rank for " << taskEligibility.first;
        EXPECT_THAT(assignments, Each(IsEmpty()));

        EXPECT_FALSE(taskAssigner_->didUserAssignGpus());
        // TODO There could be log output that we might assert upon
    }
}

INSTANTIATE_TEST_CASE_P(AllTaskEligibilitiesWork, NonbondedOnCpuTest, Combine
                            (Bool(),
                            Bool()
                            ));

} // namespace
} // namespace
} // namespace
