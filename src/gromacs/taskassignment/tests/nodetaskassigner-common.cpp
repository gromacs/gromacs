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
 * Defines helper functionality for refdata for node task-assignment tests.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "nodetaskassigner-common.h"

#include <string>
#include <vector>

#include "testutils/refdata.h"

namespace gmx
{

void gpuTaskItemChecker(test::TestReferenceChecker *checker, const GpuTask &task)
{
    if (task == GpuTask::Nonbonded)
    {
        checker->checkString("Nonbonded", "Type");
    }
    else
    {
        checker->checkString("Unknown", "Type");
    }
}

void gpuTasksItemChecker(test::TestReferenceChecker *checker, const std::vector<GpuTask> &tasks)
{
    if (tasks.empty())
    {
        checker->checkString("No GPU tasks", nullptr);
    }
    else
    {
        checker->checkSequence(tasks.begin(), tasks.end(), nullptr, gpuTaskItemChecker);
    }
}

void gpuTaskAssignmentItemChecker(test::TestReferenceChecker *checker, const GpuTaskMapping &taskMapping)
{
    gpuTaskItemChecker(checker, taskMapping.task_);
    checker->checkInteger(taskMapping.deviceId_, "Device ID");
}

void gpuTaskAssignmentsItemChecker(test::TestReferenceChecker *checker, GpuTaskAssignments::value_type &taskAssignments)
{
    checker->checkSequence(taskAssignments.begin(), taskAssignments.end(), nullptr, gpuTaskAssignmentItemChecker);
}

void gpuTaskAssignmentsOnRanksItemChecker(test::TestReferenceChecker *checker, GpuTaskAssignments &taskAssignmentsOnRanks)
{
    if (taskAssignmentsOnRanks.empty())
    {
        checker->checkString("No GPU task assignments on ranks", nullptr);
    }
    else
    {
        checker->checkSequence(taskAssignmentsOnRanks.begin(), taskAssignmentsOnRanks.end(), nullptr, gpuTaskAssignmentsItemChecker);
    }
}

} // namespace
