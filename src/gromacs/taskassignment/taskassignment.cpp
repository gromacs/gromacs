/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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

#include <algorithm>
#include <exception>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
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
std::vector<GpuTaskAssignment> buildTaskAssignment(const GpuTasksOnRanks& gpuTasksOnRanksOfThisNode,
                                                   ArrayRef<const int>    gpuIds)
{
    std::vector<GpuTaskAssignment> gpuTaskAssignmentOnRanksOfThisNode(gpuTasksOnRanksOfThisNode.size());

    // Loop over the ranks on this node, and the tasks on each
    // rank. For each task, take the next device ID from those
    // provided by the user, to build a vector of mappings of task to
    // ID, for each rank on this node. Note that if there have not
    // been any GPU tasks identified, then gpuIds can be empty.
    auto currentGpuId            = gpuIds.begin();
    auto gpuTaskAssignmentOnRank = gpuTaskAssignmentOnRanksOfThisNode.begin();
    for (const auto& gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        gpuTaskAssignmentOnRank->reserve(gpuTasksOnRank.size());
        for (const auto& gpuTaskType : gpuTasksOnRank)
        {
            GMX_RELEASE_ASSERT(currentGpuId != gpuIds.end(), "Indexing out of range for GPU tasks");
            gpuTaskAssignmentOnRank->push_back({ gpuTaskType, *currentGpuId });
            ++currentGpuId;
        }
        GMX_RELEASE_ASSERT(gpuTaskAssignmentOnRank->size() == gpuTasksOnRank.size(),
                           "Mismatch in number of GPU tasks on a rank with the number of elements "
                           "in the resulting task assignment");
        ++gpuTaskAssignmentOnRank;
    }

    return gpuTaskAssignmentOnRanksOfThisNode;
}

/*! \brief Return whether a GPU device is shared between any ranks.
 *
 * Sharing GPUs among multiple ranks is possible via either user or
 * automated selection. */
bool isAnyGpuSharedBetweenRanks(ArrayRef<const GpuTaskAssignment> gpuTaskAssignments)
{
    // Loop over all ranks i, looking on all higher ranks j whether
    // any tasks on them share GPU device IDs.
    //
    // TODO Should this functionality also consider whether tasks on
    // the same rank are sharing a device?
    for (size_t i = 0; i < gpuTaskAssignments.size(); ++i)
    {
        for (const auto& taskOnRankI : gpuTaskAssignments[i])
        {
            for (size_t j = i + 1; j < gpuTaskAssignments.size(); ++j)
            {
                for (const auto& taskOnRankJ : gpuTaskAssignments[j])
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

} // namespace

void GpuTaskAssignments::logPerformanceHints(const MDLogger& mdlog, size_t numCompatibleGpusOnThisNode)
{
    if (numCompatibleGpusOnThisNode > numGpuTasksOnThisNode_)
    {
        /* TODO In principle, this warning could be warranted only on
         * some nodes, but we lack the infrastructure to do a good job
         * of reporting that. */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: You assigned the GPU tasks on a node such that some GPUs "
                        "available on that node are unused, which might not be optimal.");
    }

    if (isAnyGpuSharedBetweenRanks(assignmentForAllRanksOnThisNode_))
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: You assigned the same GPU ID(s) to multiple ranks, which is a good "
                        "idea if you have measured the performance of alternatives.");
    }
}

namespace
{

//! Counts all the GPU tasks on this node.
size_t countGpuTasksOnThisNode(const GpuTasksOnRanks& gpuTasksOnRanksOfThisNode)
{
    size_t numGpuTasksOnThisNode = 0;
    for (const auto& gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        numGpuTasksOnThisNode += gpuTasksOnRank.size();
    }
    return numGpuTasksOnThisNode;
}

/*! \brief Return on each rank the total count over all ranks of all
 * simulations. */
int countOverAllRanks(const t_commrec* cr, const gmx_multisim_t* ms, const int countOnThisRank)
{
    int countOverAllRanksValue = countOnThisRank;
    if (PAR(cr))
    {
        // Count over the ranks of this simulation.
        gmx_sumi(1, &countOverAllRanksValue, cr);
    }
    if (isMultiSim(ms))
    {
        // Count over the ranks of all simulations.
        gmx_sumi_sim(1, &countOverAllRanksValue, ms);
        if (PAR(cr))
        {
            // Propagate the information from other simulations back
            // to non-master ranks so they can all agree on future
            // behavior.
            gmx_bcast(sizeof(decltype(countOverAllRanksValue)), &countOverAllRanksValue, cr);
        }
    }
    return countOverAllRanksValue;
}

} // namespace

GpuTaskAssignmentsBuilder::GpuTaskAssignmentsBuilder() = default;

GpuTaskAssignments GpuTaskAssignmentsBuilder::build(const std::vector<int>& gpuIdsToUse,
                                                    const std::vector<int>& userGpuTaskAssignment,
                                                    const gmx_hw_info_t&    hardwareInfo,
                                                    const t_commrec*        cr,
                                                    const gmx_multisim_t*   ms,
                                                    const PhysicalNodeCommunicator& physicalNodeComm,
                                                    const TaskTarget                nonbondedTarget,
                                                    const TaskTarget                pmeTarget,
                                                    const TaskTarget                bondedTarget,
                                                    const TaskTarget                updateTarget,
                                                    const bool useGpuForNonbonded,
                                                    const bool useGpuForPme,
                                                    bool       rankHasPpTask,
                                                    bool       rankHasPmeTask)
{
    size_t               numRanksOnThisNode = physicalNodeComm.size_;
    std::vector<GpuTask> gpuTasksOnThisRank = findGpuTasksOnThisRank(
            !gpuIdsToUse.empty(), nonbondedTarget, pmeTarget, bondedTarget, updateTarget,
            useGpuForNonbonded, useGpuForPme, rankHasPpTask, rankHasPmeTask);
    /* Communicate among ranks on this node to find each task that can
     * be executed on a GPU, on each rank. */
    auto gpuTasksOnRanksOfThisNode = findAllGpuTasksOnThisNode(gpuTasksOnThisRank, physicalNodeComm);
    size_t numGpuTasksOnThisNode   = countGpuTasksOnThisNode(gpuTasksOnRanksOfThisNode);

    std::exception_ptr             exceptionPtr;
    std::vector<GpuTaskAssignment> taskAssignmentOnRanksOfThisNode;
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

            // enforce the single device/rank restriction
            if (numRanksOnThisNode == 1 && !compatibleGpusToUse.empty())
            {
                compatibleGpusToUse = compatibleGpusToUse.subArray(0, 1);
            }

            // When doing automated assignment of GPU tasks to GPU
            // IDs, even if we have more than one kind of GPU task, we
            // do a simple round-robin assignment. That's not ideal,
            // but we don't have any way to do a better job reliably.
            generatedGpuIds = makeGpuIds(compatibleGpusToUse, numGpuTasksOnThisNode);

            if ((numGpuTasksOnThisNode > gpuIdsToUse.size())
                && (numGpuTasksOnThisNode % gpuIdsToUse.size() != 0))
            {
                // TODO Decorating the message with hostname should be
                // the job of an error-reporting module.
                char host[STRLEN];
                gmx_gethostname(host, STRLEN);

                GMX_THROW(InconsistentInputError(formatString(
                        "There were %zu GPU tasks found on node %s, but %zu GPUs were "
                        "available. If the GPUs are equivalent, then it is usually best "
                        "to have a number of tasks that is a multiple of the number of GPUs. "
                        "You should reconsider your GPU task assignment, "
                        "number of ranks, or your use of the -nb, -pme, and -npme options, "
                        "perhaps after measuring the performance you can get.",
                        numGpuTasksOnThisNode, host, gpuIdsToUse.size())));
            }
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

                GMX_THROW(InconsistentInputError(formatString(
                        "There were %zu GPU tasks assigned on node %s, but %zu GPU tasks were "
                        "identified, and these must match. Reconsider your GPU task assignment, "
                        "number of ranks, or your use of the -nb, -pme, and -npme options.",
                        userGpuTaskAssignment.size(), host, numGpuTasksOnThisNode)));
            }
            // Did the user choose compatible GPUs?
            checkUserGpuIds(hardwareInfo.gpu_info, gpuIdsToUse, userGpuTaskAssignment);

            gpuIdsForTaskAssignment = userGpuTaskAssignment;
        }
        taskAssignmentOnRanksOfThisNode =
                buildTaskAssignment(gpuTasksOnRanksOfThisNode, gpuIdsForTaskAssignment);
    }
    catch (...)
    {
        exceptionPtr = std::current_exception();
    }
    int countOfExceptionsOnThisRank   = int(bool(exceptionPtr));
    int countOfExceptionsOverAllRanks = countOverAllRanks(cr, ms, countOfExceptionsOnThisRank);

    // Avoid all ranks spamming the error stream
    //
    // TODO improve this so that unique errors on different ranks
    // are all reported on one rank.
    if (countOfExceptionsOnThisRank > 0 && physicalNodeComm.rank_ == 0)
    {
        try
        {
            if (exceptionPtr)
            {
                std::rethrow_exception(exceptionPtr);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    // TODO This implements a global barrier so that MPI runtimes can
    // organize an orderly shutdown if one of the ranks has had to
    // issue a fatal error above. When we have MPI-aware error
    // handling and reporting, this should be improved (perhaps
    // centralized there).
    simulationBarrier(cr);
    multiSimBarrier(ms);
    simulationBarrier(cr);
    if (countOfExceptionsOverAllRanks > 0)
    {
        gmx_fatal(FARGS,
                  "Exiting because task assignment failed. If there is no descriptive error "
                  "message in the terminal output, please report this failure as a bug.");
    }

    // TODO There is no check that mdrun -nb gpu or -pme gpu or
    // -gpu_id is actually being implemented such that nonbonded tasks
    // are being run on compatible GPUs, on all applicable ranks. That
    // would require communication.

    GpuTaskAssignments gpuTaskAssignments(hardwareInfo);
    gpuTaskAssignments.assignmentForAllRanksOnThisNode_ = taskAssignmentOnRanksOfThisNode;
    gpuTaskAssignments.indexOfThisRank_                 = physicalNodeComm.rank_;
    gpuTaskAssignments.numGpuTasksOnThisNode_           = numGpuTasksOnThisNode;
    gpuTaskAssignments.numRanksOnThisNode_              = numRanksOnThisNode;
    return gpuTaskAssignments;
}

GpuTaskAssignments::GpuTaskAssignments(const gmx_hw_info_t& hardwareInfo) :
    hardwareInfo_(hardwareInfo)
{
}

void GpuTaskAssignments::reportGpuUsage(const MDLogger& mdlog,
                                        bool            printHostName,
                                        bool            useGpuForBonded,
                                        PmeRunMode      pmeRunMode)
{
    gmx::reportGpuUsage(mdlog, assignmentForAllRanksOnThisNode_, numGpuTasksOnThisNode_,
                        numRanksOnThisNode_, printHostName, useGpuForBonded, pmeRunMode);
}

gmx_device_info_t* GpuTaskAssignments::initNonbondedDevice(const t_commrec* cr) const
{
    gmx_device_info_t*       deviceInfo        = nullptr;
    const GpuTaskAssignment& gpuTaskAssignment = assignmentForAllRanksOnThisNode_[indexOfThisRank_];

    // This works because only one task of each type per rank is currently permitted.
    auto nbGpuTaskMapping = std::find_if(gpuTaskAssignment.begin(), gpuTaskAssignment.end(),
                                         hasTaskType<GpuTask::Nonbonded>);
    if (nbGpuTaskMapping != gpuTaskAssignment.end())
    {
        int deviceId = nbGpuTaskMapping->deviceId_;
        deviceInfo   = getDeviceInfo(hardwareInfo_.gpu_info, deviceId);
        init_gpu(deviceInfo);

        // TODO Setting up this sharing should probably part of
        // init_domain_decomposition after further refactoring.
        if (DOMAINDECOMP(cr))
        {
            /* When we share GPUs over ranks, we need to know this for the DLB */
            dd_setup_dlb_resource_sharing(cr, deviceId);
        }
    }
    return deviceInfo;
}

gmx_device_info_t* GpuTaskAssignments::initPmeDevice() const
{
    gmx_device_info_t*       deviceInfo        = nullptr;
    const GpuTaskAssignment& gpuTaskAssignment = assignmentForAllRanksOnThisNode_[indexOfThisRank_];

    // This works because only one task of each type is currently permitted.
    auto       pmeGpuTaskMapping = std::find_if(gpuTaskAssignment.begin(), gpuTaskAssignment.end(),
                                          hasTaskType<GpuTask::Pme>);
    const bool thisRankHasPmeGpuTask = (pmeGpuTaskMapping != gpuTaskAssignment.end());
    if (thisRankHasPmeGpuTask)
    {
        deviceInfo = getDeviceInfo(hardwareInfo_.gpu_info, pmeGpuTaskMapping->deviceId_);
        init_gpu(deviceInfo);
    }
    return deviceInfo;
}

bool GpuTaskAssignments::thisRankHasPmeGpuTask() const
{
    const GpuTaskAssignment& gpuTaskAssignment = assignmentForAllRanksOnThisNode_[indexOfThisRank_];

    auto       pmeGpuTaskMapping = std::find_if(gpuTaskAssignment.begin(), gpuTaskAssignment.end(),
                                          hasTaskType<GpuTask::Pme>);
    const bool thisRankHasPmeGpuTask = (pmeGpuTaskMapping != gpuTaskAssignment.end());

    return thisRankHasPmeGpuTask;
}

bool GpuTaskAssignments::thisRankHasAnyGpuTask() const
{
    const GpuTaskAssignment& gpuTaskAssignment = assignmentForAllRanksOnThisNode_[indexOfThisRank_];

    const bool thisRankHasAnyGpuTask = !gpuTaskAssignment.empty();
    return thisRankHasAnyGpuTask;
}

} // namespace gmx
