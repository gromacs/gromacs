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
 * Defines helper and factory functionality for task assignment.
 *
 * Note that the GPU ID assignment could potentially match many
 * different kinds of simulation setups, including ranks from multiple
 * simulations, ranks from the same simulation, and/or ranks with duty
 * only for particular tasks (e.g. PME-only ranks). Which GPU ID
 * assignments are valid will naturally depend on the other run-time
 * options given to mdrun, and the current capabilities of the
 * implementation.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "taskassignment.h"

#include "config.h"

#include <string>
#include <vector>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#include "findallgputasks.h"
#include "reportgpuusage.h"

namespace gmx
{

namespace
{

/*! \brief Build data structure of types of GPU tasks on a rank,
 * together with the mapped GPU device IDs, for all GPU tasks on all
 * the ranks of this node.
 *
 * \param[in]   gpuTasksOnRanksOfThisNode  For each rank on this node, the set of tasks
 *                                         that are eligible to run on GPUs.
 * \param[in]   gpuIds                     The user-supplied GPU IDs.
 */
static GpuTaskAssignments
buildTaskAssignment(const GpuTasksOnRanks  &gpuTasksOnRanksOfThisNode,
                    ArrayRef<const int>     gpuIds)
{
    GpuTaskAssignments gpuTaskAssignmentOnRanksOfThisNode(gpuTasksOnRanksOfThisNode.size());

    // Loop over the ranks on this node, and the tasks on each
    // rank. For each task, take the next device ID from those
    // provided by the user, to build a vector of mappings of task to
    // ID, for each rank on this node. Note that if there have not
    // been any GPU tasks identified, then gpuIds can be empty.
    auto               currentGpuId            = gpuIds.begin();
    auto               gpuTaskAssignmentOnRank = gpuTaskAssignmentOnRanksOfThisNode.begin();
    for (const auto &gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        gpuTaskAssignmentOnRank->reserve(gpuTasksOnRank.size());
        for (const auto &gpuTaskType : gpuTasksOnRank)
        {
            GMX_RELEASE_ASSERT(currentGpuId != gpuIds.end(), "Indexing out of range for GPU tasks");
            gpuTaskAssignmentOnRank->push_back({gpuTaskType, *currentGpuId});
            ++currentGpuId;
        }
        GMX_RELEASE_ASSERT(gpuTaskAssignmentOnRank->size() == gpuTasksOnRank.size(),
                           "Mismatch in number of GPU tasks on a rank with the number of elements in the resulting task assignment");
        ++gpuTaskAssignmentOnRank;
    }

    return gpuTaskAssignmentOnRanksOfThisNode;
}

/*! \brief Return whether a GPU device is shared between any ranks.
 *
 * Sharing GPUs among multiple ranks is possible via either user or
 * automated selection. */
static bool isAnyGpuSharedBetweenRanks(const GpuTaskAssignments &gpuTaskAssignments)
{
    // Loop over all ranks i, looking on all higher ranks j whether
    // any tasks on them share GPU device IDs.
    //
    // TODO Should this functionality also consider whether tasks on
    // the same rank are sharing a device?
    for (size_t i = 0; i < gpuTaskAssignments.size(); ++i)
    {
        for (const auto &taskOnRankI : gpuTaskAssignments[i])
        {
            for (size_t j = i+1; j < gpuTaskAssignments.size(); ++j)
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

//! Logs to \c mdlog information that may help a user learn how to let mdrun make a task assignment that runs faster.
void logPerformanceHints(const MDLogger           &mdlog,
                         size_t                    numCompatibleGpus,
                         size_t                    numGpuTasksOnThisNode,
                         const GpuTaskAssignments &gpuTaskAssignments)
{
    if (numCompatibleGpus > numGpuTasksOnThisNode)
    {
        /* TODO In principle, this warning could be warranted only on
         * some nodes, but we lack the infrastructure to do a good job
         * of reporting that. */
        GMX_LOG(mdlog.warning).asParagraph().
            appendText("NOTE: You assigned the GPU tasks on a node such that some GPUs "
                       "available on that node are unused, which might not be optimal.");
    }

    if (isAnyGpuSharedBetweenRanks(gpuTaskAssignments))
    {
        GMX_LOG(mdlog.warning).asParagraph().
            appendText("NOTE: You assigned the same GPU ID(s) to multiple ranks, which is a good idea if you have measured the performance of alternatives.");
    }
}

//! Counts all the GPU tasks on this node.
size_t countGpuTasksOnThisNode(const GpuTasksOnRanks &gpuTasksOnRanksOfThisNode)
{
    size_t numGpuTasksOnThisNode = 0;
    for (const auto &gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        numGpuTasksOnThisNode += gpuTasksOnRank.size();
    }
    return numGpuTasksOnThisNode;
}

//! Finds whether there is any task of \c queryTask in the tasks on the ranks of this node.
bool hasAnyTaskOfTypeOnThisNode(const GpuTasksOnRanks &gpuTasksOnRanksOfThisNode,
                                const GpuTask          queryTask)
{
    for (const auto &gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        for (const auto &gpuTask : gpuTasksOnRank)
        {
            if (queryTask == gpuTask)
            {
                return true;
            }
        }
    }
    return false;
}

}   // namespace

GpuTaskAssignments::value_type
runTaskAssignment(const std::vector<int>     &gpuIdsToUse,
                  const std::vector<int>     &userGpuTaskAssignment,
                  const gmx_hw_info_t        &hardwareInfo,
                  const MDLogger             &mdlog,
                  const t_commrec            *cr,
                  const std::vector<GpuTask> &gpuTasksOnThisRank)
{
    /* Communicate among ranks on this node to find each task that can
     * be executed on a GPU, on each rank. */
    auto gpuTasksOnRanksOfThisNode = findAllGpuTasksOnThisNode(gpuTasksOnThisRank,
                                                               cr->nrank_intranode,
                                                               cr->mpi_comm_physicalnode);
    auto               numGpuTasksOnThisNode = countGpuTasksOnThisNode(gpuTasksOnRanksOfThisNode);

    GpuTaskAssignments taskAssignmentOnRanksOfThisNode;
    try
    {
        // Use the GPU IDs from the user if they supplied
        // them. Otherwise, choose from the compatible GPUs.
        //
        // GPU ID assignment strings, if provided, cover all the ranks
        // on a node. If nodes or the process placement on them are
        // heterogeneous, then the GMX_GPU_ID environment variable
        // must be set by a user who also wishes to direct GPU ID
        // assignment.  Thus this implementation of task assignment
        // can assume it has a GPU ID assignment appropriate for the
        // node upon which its process is running.
        //
        // Valid GPU ID assignments are `an ordered set of digits that
        // identify GPU device IDs (e.g. as understood by the GPU
        // runtime, and subject to environment modification such as
        // with CUDA_VISIBLE_DEVICES) that will be used for the
        // GPU-suitable tasks on all of the ranks of that node.
        ArrayRef<const int> gpuIdsForTaskAssignment;
        std::vector<int>    generatedGpuIds;
        if (userGpuTaskAssignment.empty())
        {
            ArrayRef<const int> compatibleGpusToUse = gpuIdsToUse;
            if (hasAnyTaskOfTypeOnThisNode(gpuTasksOnRanksOfThisNode, GpuTask::Pme))
            {
                // PP and PME tasks must run on the same device, so
                // restrict the assignment to the first device. If
                // there aren't any, then that error is handled later.
                if (!compatibleGpusToUse.empty())
                {
                    compatibleGpusToUse = compatibleGpusToUse.subArray(0, 1);
                }
            }
            generatedGpuIds         = makeGpuIds(compatibleGpusToUse, numGpuTasksOnThisNode);
            gpuIdsForTaskAssignment = generatedGpuIds;
        }
        else
        {
            if (numGpuTasksOnThisNode != userGpuTaskAssignment.size())
            {
                // TODO Decorating the message with hostname should be
                // the job of an error-reporting module.
                char host[STRLEN];
                gmx_gethostname(host, STRLEN);

                GMX_THROW(InconsistentInputError
                              (formatString("There were %zu GPU tasks assigned on node %s, but %zu GPU tasks were "
                                            "identified, and these must match. Reconsider your GPU task assignment, "
                                            "number of ranks, or your use of the -nb, -pme, and -npme options.", userGpuTaskAssignment.size(),
                                            host, numGpuTasksOnThisNode)));
            }
            // Did the user choose compatible GPUs?
            checkUserGpuIds(hardwareInfo.gpu_info, gpuIdsToUse, userGpuTaskAssignment);

            gpuIdsForTaskAssignment = userGpuTaskAssignment;
        }
        taskAssignmentOnRanksOfThisNode =
            buildTaskAssignment(gpuTasksOnRanksOfThisNode, gpuIdsForTaskAssignment);

    }
    catch (const std::exception &ex)
    {
        // TODO This implementation is quite similar to that of
        // processExceptionAsFatalError (which implements
        // GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR), but it is unclear
        // how we should involve MPI in the implementation of error
        // handling.
        if (cr->rank_intranode == 0)
        {
            printFatalErrorMessage(stderr, ex);
        }

        if (PAR(cr))
        {
#if GMX_MPI
            MPI_Barrier(cr->mpi_comm_mysim);
#endif
        }
        if (MULTISIM(cr))
        {
#if GMX_MPI
            MPI_Barrier(cr->ms->mpi_comm_masters);
#endif
        }

        gmx_exit_on_fatal_error(ExitType_Abort, 1);
    }

    reportGpuUsage(mdlog, !userGpuTaskAssignment.empty(), taskAssignmentOnRanksOfThisNode,
                   numGpuTasksOnThisNode, cr->nrank_intranode, cr->nnodes > 1);

    // If the user chose a task assignment, give them some hints where appropriate.
    if (!userGpuTaskAssignment.empty())
    {
        logPerformanceHints(mdlog, gpuIdsToUse.size(),
                            numGpuTasksOnThisNode,
                            taskAssignmentOnRanksOfThisNode);
    }

    return taskAssignmentOnRanksOfThisNode[cr->rank_intranode];

    // TODO There is no check that mdrun -nb gpu or -pme gpu or
    // -gpu_id is actually being implemented such that nonbonded tasks
    // are being run on compatible GPUs, on all applicable ranks. That
    // would require communication.
}

} // namespace
