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
 * \brief Defines routine for reporting GPU usage.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "reportgpuusage.h"

#include <set>
#include <string>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

namespace gmx
{

namespace
{

/*! \brief Count and return the number of unique GPUs (per node) selected.
 *
 * As sharing GPUs among multiple ranks is possible, the number of
 * GPUs used (per node) can be different from the number of GPU IDs
 * used.
 */
static size_t countUniqueGpuIdsUsed(const GpuTaskAssignments &gpuTaskAssignmentOnRanksOfThisNode)
{
    std::set<int> uniqueIds;
    for (const auto &assignmentsOnRank : gpuTaskAssignmentOnRanksOfThisNode)
    {
        for (const auto &assignmentOfTask : assignmentsOnRank)
        {
            uniqueIds.insert(assignmentOfTask.deviceId_);
        }
    }
    return uniqueIds.size();
}

}   // namespace

void
reportGpuUsage(const MDLogger                &mdlog,
               bool                           userSetGpuIds,
               const GpuTaskAssignments      &gpuTaskAssignmentOnRanksOfThisNode,
               size_t                         numGpuTasksOnThisNode,
               size_t                         numRanks,
               bool                           bPrintHostName)
{
    size_t numGpusInUse = countUniqueGpuIdsUsed(gpuTaskAssignmentOnRanksOfThisNode);
    if (numGpusInUse == 0)
    {
        return;
    }

    std::string output;
    {
        std::string gpuIdsString;
        const char *currentSeparator = "";
        const char *separator        = ",";
        for (const auto &assignmentsOnRank : gpuTaskAssignmentOnRanksOfThisNode)
        {
            if (assignmentsOnRank.empty())
            {
                gpuIdsString    += currentSeparator;
                gpuIdsString    += "none";
                currentSeparator = separator;
            }
            else
            {
                for (const auto &assignmentOnRank : assignmentsOnRank)
                {
                    const char *rankType = (assignmentOnRank.task_ == GpuTask::Nonbonded ? "PP" : "PME");
                    gpuIdsString    += currentSeparator;
                    gpuIdsString    += formatString("%s:%d", rankType, assignmentOnRank.deviceId_);
                    currentSeparator = separator;
                }
            }
        }
        bool        bPluralGpus  = numGpusInUse > 1;

        if (bPrintHostName)
        {
            char host[STRLEN];
            gmx_gethostname(host, STRLEN);
            output += gmx::formatString("On host %s ", host);
        }
        output += gmx::formatString("%zu GPU%s %sselected for this run.\n"
                                    "Mapping of GPU IDs to the %d GPU task%s in the %d rank%s on this node:\n  %s\n",
                                    numGpusInUse, bPluralGpus ? "s" : "",
                                    userSetGpuIds ? "user-" : "auto-",
                                    numGpuTasksOnThisNode,
                                    (numGpuTasksOnThisNode > 1) ? "s" : "",
                                    numRanks,
                                    (numRanks > 1) ? "s" : "",
                                    gpuIdsString.c_str());
    }

    /* NOTE: this print is only for and on one physical node */
    GMX_LOG(mdlog.warning).appendText(output);
}

} // namespace
