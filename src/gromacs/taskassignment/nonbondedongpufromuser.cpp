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
 * Defines gmx::INodeTaskAssigner interface, its concrete classes, and factory function.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "nonbondedongpufromuser.h"

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

struct gmx_gpu_info_t;

namespace gmx
{

NonbondedOnGpuFromUser::NonbondedOnGpuFromUser(const std::string   &gpuIdString,
                                               const gmx_hw_info_t &hardwareInfo)
    : gpuIdString_(gpuIdString), hardwareInfo_(hardwareInfo) {}

void
NonbondedOnGpuFromUser::decideWhetherToUseGpus(const bool usingVerletScheme,
                                               const bool gpuAccelerationIsUseful)
{
    if (!usingVerletScheme)
    {
        GMX_THROW(InconsistentInputError
                      ("GPU acceleration required, but can't be used without cutoff-scheme=Verlet."));
    }
    if (!gpuAccelerationIsUseful)
    {
        // Details about why it isn't useful should have been issued
        // as a warning by the code that realized that, because it's
        // only a fatal error for this INodeTaskAssigner.
        GMX_THROW(InconsistentInputError
                      ("GPU acceleration required, but that is not supported with the given input "
                      "settings. Either do not request GPUs, or change your inputs to something supported."));
    }
}

bool NonbondedOnGpuFromUser::areNonbondedOnGpu() const
{
    return true;
}

namespace
{

/*! \brief Check that all user-selected GPUs are compatible.
 *
 * Given the \c userGpuTaskAssignment and \c hardwareInfo, throw if
 * any selected GPUs is not compatible.
 *
 * The error is given with a suitable descriptive message, which will
 * have context if this check is done after the hardware detection
 * results have been reported to the user. However, note that only the
 * GPUs detected on the master rank are reported, because of the
 * existing limitations of that reporting.
 *
 * \todo Note that the selected GPUs can be different on each rank,
 * and the IDs of compatible GPUs can be different on each node, so
 * this routine ought to do communication to determine whether all
 * ranks are able to proceed. Currently this relies on the MPI runtime
 * to kill the other processes because GROMACS lacks the appropriate
 * infrastructure to do a good job of coordinating error messages and
 * behaviour across MPMD ranks and multiple simulations.
 *
 * \param[in]   hardwareInfo           The detected hardware
 * \param[in]   userGpuTaskAssignment  The GPU selection from the user.
 *
 * \throws  std::bad_alloc          If out of memory
 *          InconsistentInputError  If the assigned GPUs are not valid
 */
static void throwUnlessUserGpuTaskAssignmentIsValid(const gmx_hw_info_t    &hardwareInfo,
                                                    const std::vector<int> &userGpuTaskAssignment)
{
    bool        foundIncompatibleGpuIds = false;
    std::string message
        = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    const auto &compatibleGpus = hardwareInfo.compatibleGpus;
    for (const auto &gpuId : userGpuTaskAssignment)
    {
        if (std::find(compatibleGpus.begin(), compatibleGpus.end(), gpuId) == compatibleGpus.end())
        {
            foundIncompatibleGpuIds = true;
            message                += gmx::formatString("    GPU #%d: %s\n",
                                                        gpuId,
                                                        getGpuCompatibilityDescription(hardwareInfo.gpu_info, gpuId));
        }
    }

    if (foundIncompatibleGpuIds)
    {
        GMX_THROW(InconsistentInputError(message));
    }
}

}   // namespace

GpuTaskAssignments
NonbondedOnGpuFromUser::assignGpuTasksToDeviceIds(const GpuTasksOnRanks &gpuTasksOnRanks) const
{
    auto numGpuTasksOnThisNode = countGpuTasksOnThisNode(gpuTasksOnRanks);

    if (numGpuTasksOnThisNode == 0)
    {
        // No GPU tasks found, so return an empty assignment. This
        // could be a PME-only rank, which ignores -gpu_id. Even if
        // the user required the use of GPUs for nonbondeds, the check
        // that some node has ranks that will use GPUs for nonbondeds
        // must take place elsewhere.
        return GpuTaskAssignments(gpuTasksOnRanks.size());
    }

    char host[STRLEN];
    gmx_gethostname(host, STRLEN);

    if (numGpuTasksOnThisNode != gpuIdString_.size())
    {
        auto message = formatString("You selected %lu GPU IDs for the GPU tasks on node %s "
                                    "but %lu GPU tasks were found on the ranks of this node",
                                    gpuIdString_.size(), host, numGpuTasksOnThisNode);
        GMX_THROW(InconsistentInputError(message));
    }
    /* Validate the user-selected assignment of short-ranged tasks
     * on PP ranks (from mdrun -gpu_id), now that it is known from
     * the tasks on each rank whether the available GPU tasks
     * match the GPU IDs the user selected. */
    auto userGpuTaskAssignment = parseGpuTaskAssignment(gpuIdString_);
    throwUnlessUserGpuTaskAssignmentIsValid(hardwareInfo_, userGpuTaskAssignment);

    return buildTaskAssignment(gpuTasksOnRanks, userGpuTaskAssignment);
}

bool
NonbondedOnGpuFromUser::didUserAssignGpus() const
{
    return true;
}

namespace
{

/*! \brief Return whether a GPU device is shared between any ranks.
 *
 * Sharing GPUs among multiple ranks is possible via either user or
 * automated selection. */
static bool isAnyGpuSharedBetweenRanks(const GpuTaskAssignments &gpuTaskAssignments)
{
    /* Loop over all ranks i, looking on all higher ranks j whether
       any tasks on them share GPU device IDs. */
    for (size_t i = 0; i != gpuTaskAssignments.size(); ++i)
    {
        for (const auto &taskOnRankI : gpuTaskAssignments[i])
        {
            for (size_t j = 0; j != gpuTaskAssignments.size(); ++j)
            {
                for (const auto &taskOnRankJ : gpuTaskAssignments[j])
                {
                    if (taskOnRankI.deviceId_ == taskOnRankJ.deviceId_)
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

}   // namespace

void
NonbondedOnGpuFromUser::logPerformanceHints(const MDLogger           &mdlog,
                                            const GpuTaskAssignments &gpuTaskAssignments) const
{
    auto numGpuTasksOnThisNode = countGpuTasksOnThisNode(gpuTaskAssignments);

    if (hardwareInfo_.compatibleGpus.size() > numGpuTasksOnThisNode)
    {
        /* TODO In principle, this warning could be warranted only on
         * some nodes, but we lack the infrastructure to do a good job
         * of reporting that. */
        GMX_LOG(mdlog.warning).asParagraph().
            appendText("NOTE: You assigned the GPU tasks on a node such that some GPUs"
                       "available on that node are unused, which might not be optimal.");
    }

    if (isAnyGpuSharedBetweenRanks(gpuTaskAssignments))
    {
        GMX_LOG(mdlog.warning).asParagraph().
            appendText("NOTE: You assigned the same GPU ID(s) to multiple ranks, which is a good idea if you have measured the performance of alternatives.");
    }

}

} // namespace
