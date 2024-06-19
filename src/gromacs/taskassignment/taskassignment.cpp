/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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

#include "gromacs/taskassignment/taskassignment.h"

#include "config.h"

#include <cstddef>

#include <algorithm>
#include <exception>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/utility/arrayref.h"
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

enum class PmeRunMode;
struct DeviceInformation;

namespace gmx
{
enum class TaskTarget;

namespace
{

/*! \brief Build the GPU task assignment for the ranks of this node.
 *
 * \param[in]   gpuTasksOnRanksOfThisNode  For each rank on this node, the set of tasks
 *                                         that are eligible to run on GPUs.
 * \param[in]   gpuIds                     The GPU IDs for the tasks on this node, supplied
 *                                         either by the user or the automatic assignment.
 * \return      A vector with elements for each rank on this node that
 *              describes the GPU tasks and the assigned device ID.
 *
 * \throws InvalidInputError  when the user GPU assignment requests multiple devices on a rank
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
        if (gpuTasksOnRank.empty())
        {
            continue;
        }
        gpuTaskAssignmentOnRank->reserve(gpuTasksOnRank.size());
        // Keep a copy of the first GPU ID for this rank so we can
        // check that it is the only GPU ID used on this rank.
        const int gpuIdOnRank = *currentGpuId;
        for (const auto& gpuTaskType : gpuTasksOnRank)
        {
            GMX_RELEASE_ASSERT(currentGpuId != gpuIds.end(), "Indexing out of range for GPU tasks");
            gpuTaskAssignmentOnRank->push_back({ gpuTaskType, gpuIdOnRank });
            if (*currentGpuId != gpuIdOnRank)
            {
                const char* message =
                        "The GPU task assignment requested mdrun to use more than one GPU device "
                        "on a rank, which is not supported. Request only one GPU device per rank.";
                GMX_THROW(InvalidInputError(message));
            }
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

void GpuTaskAssignments::logPerformanceHints(const MDLogger& mdlog, size_t numAvailableDevicesOnThisNode)
{
    if (numAvailableDevicesOnThisNode > numGpuTasksOnThisNode_)
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
int countOverAllRanks(MPI_Comm comm, int countOnThisRank)
{
    int sum;
#if GMX_MPI
    int numRanks;
    MPI_Comm_size(comm, &numRanks);
    if (numRanks > 1)
    {
        MPI_Allreduce(&countOnThisRank, &sum, 1, MPI_INT, MPI_SUM, comm);
    }
    else
#else
    GMX_UNUSED_VALUE(comm);
#endif
    {
        sum = countOnThisRank;
    }

    return sum;
}

/*! \brief Barrier over all rank in \p comm */
void barrierOverAllRanks(MPI_Comm comm)
{
#if GMX_MPI
    int numRanks;
    MPI_Comm_size(comm, &numRanks);
    if (numRanks > 1)
    {
        MPI_Barrier(comm);
    }
#else
    GMX_UNUSED_VALUE(comm);
#endif
}

} // namespace

GpuTaskAssignmentsBuilder::GpuTaskAssignmentsBuilder() = default;

GpuTaskAssignments GpuTaskAssignmentsBuilder::build(const gmx::ArrayRef<const int> availableDevices,
                                                    const gmx::ArrayRef<const int> userGpuTaskAssignment,
                                                    const gmx_hw_info_t&           hardwareInfo,
                                                    MPI_Comm                       gromacsWorldComm,
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
    std::vector<GpuTask> gpuTasksOnThisRank = findGpuTasksOnThisRank(!availableDevices.empty(),
                                                                     nonbondedTarget,
                                                                     pmeTarget,
                                                                     bondedTarget,
                                                                     updateTarget,
                                                                     useGpuForNonbonded,
                                                                     useGpuForPme,
                                                                     rankHasPpTask,
                                                                     rankHasPmeTask);
    /* Communicate among ranks on this node to find each task that can
     * be executed on a GPU, on each rank. */
    auto gpuTasksOnRanksOfThisNode = findAllGpuTasksOnThisNode(gpuTasksOnThisRank, physicalNodeComm);
    size_t numGpuTasksOnThisNode   = countGpuTasksOnThisNode(gpuTasksOnRanksOfThisNode);

    std::exception_ptr             exceptionPtr;
    std::vector<GpuTaskAssignment> taskAssignmentOnRanksOfThisNode;
    std::vector<int>               deviceIdAssignment;
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
        std::vector<int> generatedGpuIds;
        if (userGpuTaskAssignment.empty())
        {
            ArrayRef<const int> compatibleGpusToUse = availableDevices;

            // Enforce the single device per rank restriction by ensuring
            // that there are only at most as many devices used as ranks.
            //
            // This means that with a single rank with NB and PME
            // offloaded we assign both tasks to the same GPU
            // regardless of how many GPUs are detected. Similarly,
            // with N combined PP-PME ranks (ie. with NB and PME
            // offloaded and PME decomposition active) we assign both
            // tasks on each rank to the same GPU even when more than
            // N GPUs are detected.
            if (numRanksOnThisNode < compatibleGpusToUse.size())
            {
                compatibleGpusToUse = compatibleGpusToUse.subArray(0, numRanksOnThisNode);
            }

            // When doing automated assignment of GPU tasks to GPU
            // IDs, even if we have more than one kind of GPU task, we
            // do a simple round-robin assignment. That's not ideal,
            // but we don't have any way to do a better job reliably.
            generatedGpuIds = makeGpuIds(compatibleGpusToUse, numGpuTasksOnThisNode);

            if ((numGpuTasksOnThisNode > availableDevices.size())
                && (numGpuTasksOnThisNode % availableDevices.size() != 0))
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
                        numGpuTasksOnThisNode,
                        host,
                        availableDevices.size())));
            }
            deviceIdAssignment = generatedGpuIds;
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
                        userGpuTaskAssignment.size(),
                        host,
                        numGpuTasksOnThisNode)));
            }
            // Did the user choose compatible GPUs?
            checkUserGpuIds(hardwareInfo.deviceInfoList, availableDevices, userGpuTaskAssignment);

            deviceIdAssignment = gmx::copyOf(userGpuTaskAssignment);
        }
        taskAssignmentOnRanksOfThisNode =
                buildTaskAssignment(gpuTasksOnRanksOfThisNode, deviceIdAssignment);
    }
    catch (...)
    {
        exceptionPtr = std::current_exception();
    }
    int countOfExceptionsOnThisRank = int(bool(exceptionPtr));
    int countOfExceptionsOverAllRanks = countOverAllRanks(gromacsWorldComm, countOfExceptionsOnThisRank);

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
    // TODO Global barrier so that MPI runtimes can
    // organize an orderly shutdown if one of the ranks has had to
    // issue a fatal error above. When we have MPI-aware error
    // handling and reporting, this should be improved (perhaps
    // centralized there).
    barrierOverAllRanks(gromacsWorldComm);
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
    gpuTaskAssignments.deviceIdsAssigned_               = deviceIdAssignment;
    std::sort(gpuTaskAssignments.deviceIdsAssigned_.begin(), gpuTaskAssignments.deviceIdsAssigned_.end());
    gpuTaskAssignments.deviceIdsAssigned_.erase(unique(gpuTaskAssignments.deviceIdsAssigned_.begin(),
                                                       gpuTaskAssignments.deviceIdsAssigned_.end()),
                                                gpuTaskAssignments.deviceIdsAssigned_.end());
    return gpuTaskAssignments;
}

GpuTaskAssignments::GpuTaskAssignments(const gmx_hw_info_t& hardwareInfo) :
    hardwareInfo_(hardwareInfo)
{
}

void GpuTaskAssignments::reportGpuUsage(const MDLogger&           mdlog,
                                        bool                      printHostName,
                                        PmeRunMode                pmeRunMode,
                                        const SimulationWorkload& simulationWork)
{
    gmx::reportGpuUsage(mdlog,
                        assignmentForAllRanksOnThisNode_,
                        numGpuTasksOnThisNode_,
                        numRanksOnThisNode_,
                        printHostName,
                        pmeRunMode,
                        simulationWork);
}

/*! \brief Function for whether the task of \c mapping has value \c TaskType.
 *
 * \param[in] mapping  Current GPU task mapping.
 * \returns If \c TaskType task was assigned to the \c mapping.
 */
template<GpuTask TaskType>
static bool hasTaskType(const GpuTaskMapping& mapping)
{
    return mapping.task_ == TaskType;
}

/*! \brief Function for whether the \c mapping has the GPU PME or Nonbonded task.
 *
 * \param[in] mapping  Current GPU task mapping.
 * \returns If PME on Nonbonded GPU task was assigned to this mapping.
 */
static bool hasPmeOrNonbondedTask(const GpuTaskMapping& mapping)
{
    return hasTaskType<GpuTask::Pme>(mapping) || hasTaskType<GpuTask::Nonbonded>(mapping);
}

DeviceInformation* GpuTaskAssignments::initDevice(int* deviceId) const
{
    DeviceInformation*       deviceInfo        = nullptr;
    const GpuTaskAssignment& gpuTaskAssignment = assignmentForAllRanksOnThisNode_[indexOfThisRank_];

    // This works because only one task of each type per rank is
    // currently permitted and if there are multiple tasks, they must
    // use the same device.
    auto gpuTaskMapping =
            std::find_if(gpuTaskAssignment.begin(), gpuTaskAssignment.end(), hasPmeOrNonbondedTask);

    if (gpuTaskMapping != gpuTaskAssignment.end())
    {
        *deviceId  = gpuTaskMapping->deviceId_;
        deviceInfo = hardwareInfo_.deviceInfoList[*deviceId].get();
    }
    return deviceInfo;
}

bool GpuTaskAssignments::thisRankHasPmeGpuTask() const
{
    const GpuTaskAssignment& gpuTaskAssignment = assignmentForAllRanksOnThisNode_[indexOfThisRank_];

    auto pmeGpuTaskMapping = std::find_if(
            gpuTaskAssignment.begin(), gpuTaskAssignment.end(), hasTaskType<GpuTask::Pme>);
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
