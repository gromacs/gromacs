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
 * Declares helper functionality and common test input data for node
 * task-assignment tests.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/taskassignment/inodetaskassigner.h"

#include "testutils/refdata.h"

namespace gmx
{

/**@{*/
/*! \brief Helper functions for writing refdata. */

void gpuTaskItemChecker(test::TestReferenceChecker *checker, const GpuTask &task);
void gpuTasksItemChecker(test::TestReferenceChecker *checker, const std::vector<GpuTask> &tasks);
void gpuTaskAssignmentItemChecker(test::TestReferenceChecker *checker, const GpuTaskMapping &taskMapping);
void gpuTaskAssignmentsItemChecker(test::TestReferenceChecker *checker, GpuTaskAssignments::value_type &taskAssignments);
void gpuTaskAssignmentsOnRanksItemChecker(test::TestReferenceChecker *checker, GpuTaskAssignments &taskAssignmentsOnRanks);

/**@}*/

namespace test
{
namespace
{

//! Convenience using statement.
using GpuCompatibilities = std::vector<int>;

/**@{*/
/*! \brief Values (and containers of values) used as input parameters for specifc test cases. */

std::map<std::string, GpuTasksOnRanks>     allTaskEligibilities = {
    { "One Rank No Tasks", {{}}},
    { "Two Ranks No Tasks", {{}, {}}},
    { "One Rank Each With One Nonbonded Task",
      {
          { GpuTask::Nonbonded }
      }},
    { "Two Ranks Each With One Nonbonded Task",
      {
          { GpuTask::Nonbonded },
          { GpuTask::Nonbonded }
      }},
    { "Three Ranks Each With One Nonbonded Task",
      {
          { GpuTask::Nonbonded },
          {
              GpuTask::Nonbonded
          }, {
              GpuTask::Nonbonded
          }
      }},
    { "One Rank Each With One Nonbonded Task And Rank With No Task", {{
                                                                          GpuTask::Nonbonded
                                                                      }, {}}},
    { "Two Ranks Each With One Nonbonded Task And Rank With No Task", {{
                                                                           GpuTask::Nonbonded
                                                                       }, {
                                                                           GpuTask::Nonbonded
                                                                       }, {}}},
    { "Three Ranks Each With One Nonbonded Task And Rank With No Task", {{
                                                                             GpuTask::Nonbonded
                                                                         }, {
                                                                             GpuTask::Nonbonded
                                                                         }, {
                                                                             GpuTask::Nonbonded
                                                                         }, {}}}
};

std::map<std::string, GpuCompatibilities>  allGpuCompatibilities = {
    { "No compatible devices", {}},
    { "Device 0 is compatible", {
          0
      }},
    { "Device 1 is compatible", {
          1
      }},
    { "Devices 0 and 1 are compatible", {
          0, 1
      }},
    { "Devices 0 and 2 are compatible", {
          0, 2
      }}
};

std::string
    nothing(""),
anything("anything"),
zero("0"),
one("1"),
zeroOne("01"),
zeroTwo("02");

std::vector<std::string> allUserAssignedGpuDeviceIds {
    nothing, zero, one, zeroOne, zeroTwo
};

/**@}*/

} // namespace
} // namespace
} // namespace
