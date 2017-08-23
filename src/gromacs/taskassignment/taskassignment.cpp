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
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "config.h"

#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "gromacs/compat/make_unique.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/taskassignment/inodetaskassigner.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#include "nonbondedoncpu.h"
#include "nonbondedongpuauto.h"
#include "nonbondedongpufromuser.h"

namespace gmx
{

std::vector<int> parseGpuTaskAssignment(const std::string &gpuTaskAssignment)
{
    std::vector<int> digits;
    if (gpuTaskAssignment.empty())
    {
        return digits;
    }

    /* Parse a "plain" or comma-separated GPU ID string which contains
     * a sequence of digits corresponding to GPU IDs; the order will
     * indicate the assignment of GPU tasks on this node to GPU
     * device IDs on this node. */
    try
    {
        digits = parseDigitsFromString(gpuTaskAssignment);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if (digits.empty())
    {
        gmx_fatal(FARGS, "Empty GPU ID string encountered.\n"
                  "An empty, delimiter-free, or comma-separated sequence of valid numeric IDs of available GPUs is required.\n");
    }
    return digits;
}

INodeTaskAssignerPtr
makeNodeTaskAssigner(const char          *nbOptionString,
                     const std::string   &gpuIdString,
                     EmulateGpuNonbonded  emulateGpuNonbonded,
                     const gmx_hw_info_t &hardwareInfo)
{
    bool forceNonbondedOnCpu = (strncmp(nbOptionString, "cpu", 3) == 0);

    if (forceNonbondedOnCpu)
    {
        if (!gpuIdString.empty())
        {
            GMX_THROW(InconsistentInputError("GPU IDs were specified, and short-ranged interactions were "
                                             "assigned to the CPU. Make no more than one of these choices."));
        }
        return compat::make_unique<NonbondedOnCpu>();
    }

    bool forceNonbondedOnPhysicalGpu = (strncmp(nbOptionString, "gpu", 3) == 0);
    if (emulateGpuNonbonded == EmulateGpuNonbonded::Yes)
    {
        if (forceNonbondedOnPhysicalGpu)
        {
            auto message = "Short-ranged interactions on the GPU were required, which is inconsistent "
                "with choosing emulation. Make no more than one of these choices.";
            GMX_THROW(InconsistentInputError(message));
        }
        if (!gpuIdString.empty())
        {
            auto message = "GPU IDs were specified, as was GPU emulation. Make no more than one of these choices.";
            GMX_THROW(InconsistentInputError(message));
        }
        return compat::make_unique<NonbondedOnCpu>();
    }

    // nonbonded are not forced to run the CPU (but might end up
    // there if no GPUs are detected, etc.)

    // Note that it is not (yet) an error for compatible GPUs to be
    // missing, e.g. because we have a node with only PME ranks, which
    // ignore -nb gpu and/or the short-ranged GPU IDs. So we can't
    // handle such aspects in this factory.
    if (!gpuIdString.empty())
    {
        return compat::make_unique<NonbondedOnGpuFromUser>(gpuIdString, hardwareInfo);
    }

    // TODO mdrun -nb gpu_cpu currently reaches here. What is its
    // future? If it should go away, then we can assert on the
    // contents of nbOptionString.

    // If we get here, then the user permitted automated assignment.
    return compat::make_unique<NonbondedOnGpuAuto>(forceNonbondedOnPhysicalGpu, hardwareInfo);
}

namespace
{

//! Constant used to help minimize preprocessing of code.
constexpr bool g_usingMpi = GMX_MPI;

//! Helper function to prepare to all-gather the vector of non-bonded tasks on this node.
static std::vector<int> allgather(const int &input,
                                  int        numRanks,
                                  MPI_Comm   communicator)
{
    std::vector<int> result(numRanks);
    if (g_usingMpi && numRanks > 1)
    {
        // TODO This works as an MPI_Allgather, but thread-MPI does
        // not implement that. It's only intra-node communication, and
        // happens rarely, so not worth optimizing (yet). Also
        // thread-MPI segfaults with 1 rank.
        int root = 0;
#if GMX_MPI
        MPI_Gather(&input, 1, MPI_INT,
                   result.data(), 1, MPI_INT,
                   root, communicator);
        MPI_Bcast(result.data(), result.size(), MPI_INT,
                  root, communicator);
#endif
    }
    else
    {
        result[0] = input;
    }

    return result;
}

//! Helper function to compute allgatherv displacements.
static std::vector<int> computeDisplacements(ConstArrayRef<int> extentOnEachRank,
                                             int                numRanks)
{
    std::vector<int> displacements(numRanks + 1);
    displacements[0] = 0;
    std::partial_sum(std::begin(extentOnEachRank), std::end(extentOnEachRank), std::begin(displacements) + 1);
    return displacements;
}

//! Helper function to all-gather the vector of all GPU tasks on ranks of this node.
static std::vector<GpuTask> allgatherv(ConstArrayRef<GpuTask> input,
                                       ConstArrayRef<int>     extentOnEachRank,
                                       ConstArrayRef<int>     displacementForEachRank,
                                       MPI_Comm               communicator)
{
    // Now allocate the vector and do the allgatherv
    int                  totalExtent = displacementForEachRank.back();

    std::vector<GpuTask> result(totalExtent);
    if (g_usingMpi && extentOnEachRank.size() > 1 && totalExtent > 0)
    {
        // TODO This works as an MPI_Allgatherv, but thread-MPI does
        // not implement that. It's only intra-node communication, and
        // happens rarely, so not worth optimizing (yet). Also
        // thread-MPI segfaults with 1 rank and with zero totalExtent.
        int  root = 0;
        // Calling a C API with the const int * from data() doesn't seem to compile reliably.
        auto extents       = const_cast<int *>(extentOnEachRank.data());
        auto displacements = const_cast<int *>(displacementForEachRank.data());
#if GMX_MPI
        MPI_Gatherv(input.data(), input.size(), MPI_INT,
                    result.data(), extents, displacements, MPI_INT,
                    root, communicator);
        MPI_Bcast(result.data(), result.size(), MPI_INT,
                  root, communicator);
#endif
    }
    else
    {
        for (const auto &gpuTask : input)
        {
            result.push_back(gpuTask);
        }
    }
    return result;
}

/*! \brief Returns container of all tasks on all ranks of this node
 * that are eligible for GPU execution.
 *
 * Perform all necessary communication for preparing for task
 * assignment. Separating this aspect makes it possible to unit test
 * the logic of task assignment. */
GpuTasksOnRanks
findAllGpuTasksOnThisNode(ConstArrayRef<ComputeTask> computeTasksOnThisRank,
                          bool                       nonbondedOnGpu,
                          int                        numRanksOnThisNode,
                          MPI_Comm                   communicator)
{
    /* Filter the tasks on this rank to become the vector of tasks on
     * this rank that can run on GPUs. For now, that's only non-bonded
     * tasks. */
    std::vector<GpuTask> gpuTasksOnThisRank;
    for (const auto &taskOnRank : computeTasksOnThisRank)
    {
        // Can this task run on a GPU?
        // Later, PME-on-GPU code will enable long-ranged tasks.
        if (taskOnRank == ComputeTask::Nonbonded && nonbondedOnGpu)
        {
            gpuTasksOnThisRank.push_back(GpuTask::Nonbonded);
        }
    }

    // Find out how many GPU tasks are on each rank on this node.
    auto numGpuTasksOnEachRankOfThisNode =
        allgather(gpuTasksOnThisRank.size(), numRanksOnThisNode, communicator);

    /* Collect on each rank of this node a flat vector describing all
     * GPU tasks on this node, in ascending order of rank. This
     * requires a vector allgather. The displacements indicate where
     * the GPU tasks on each rank of this node start and end within
     * the flat vector. */
    auto displacementsForEachRank = computeDisplacements(numGpuTasksOnEachRankOfThisNode, numRanksOnThisNode);
    auto flatGpuTasksOnThisNode   =
        allgatherv(gpuTasksOnThisRank, numGpuTasksOnEachRankOfThisNode,
                   displacementsForEachRank, communicator);

    /* Next, we re-use the displacements to break up the flat vector
     * of GPU tasks into something that can be indexed like
     * gpuTasks[rankIndex][taskIndex]. */
    GpuTasksOnRanks gpuTasksOnRanksOfThisNode;
    // TODO This would be nicer if we had a good abstraction for "pair
    // of iterators that point to adjacent container elements" or
    // "iterator that points to the first of a pair of valid adjacent
    // container elements, or end".
    GMX_ASSERT(displacementsForEachRank.size() > 1, "Even with one rank, there's always both a start and end displacement");
    auto currentDisplacement = displacementsForEachRank.begin();
    auto nextDisplacement = ++currentDisplacement;
    do
    {
        gpuTasksOnRanksOfThisNode.emplace_back(&flatGpuTasksOnThisNode[*currentDisplacement],
                                               &flatGpuTasksOnThisNode[*nextDisplacement]);
        currentDisplacement = nextDisplacement;
        ++nextDisplacement;
    } while (nextDisplacement != displacementsForEachRank.end());

    return gpuTasksOnRanksOfThisNode;
}

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

/*! \brief Log a report on how GPUs are being used on
 * the ranks of the physical node of rank 0 of the simulation.
 *
 * \todo It could be useful to report also whether any nodes differed,
 * and in what way.
 *
 * \param[out] mdlog                               Logging object.
 * \param[in]  userSetGpuIds                       Whether the user selected the GPU ids
 * \param[in]  gpuTaskAssignmentOnRanksOfThisNode  The selected GPU IDs.
 * \param[in]  numPpRanks                          Number of PP ranks on this node
 * \param[in]  bPrintHostName                      Print the hostname in the usage information
 *
 * \throws                        std::bad_alloc if out of memory */
static void
reportGpuUsage(const MDLogger                &mdlog,
               bool                           userSetGpuIds,
               const GpuTaskAssignments      &gpuTaskAssignmentOnRanksOfThisNode,
               size_t                         numPpRanks,
               bool                           bPrintHostName)
{
    if (gpuTaskAssignmentOnRanksOfThisNode.empty())
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
            // TODO Enhance this when there might be PP and PME tasks on a rank
            for (const auto &assignmentOnRank : assignmentsOnRank)
            {
                gpuIdsString    += currentSeparator;
                gpuIdsString    += formatString("%d", assignmentOnRank.deviceId_);
                currentSeparator = separator;
            }
        }
        size_t      numGpusInUse = countUniqueGpuIdsUsed(gpuTaskAssignmentOnRanksOfThisNode);
        bool        bPluralGpus  = numGpusInUse > 1;

        if (bPrintHostName)
        {
            char host[STRLEN];
            gmx_gethostname(host, STRLEN);
            output += gmx::formatString("On host %s ", host);
        }
        output += gmx::formatString("%zu GPU%s %sselected for this run.\n"
                                    "Mapping of GPU ID%s to the %d rank%s in this node: %s\n",
                                    numGpusInUse, bPluralGpus ? "s" : "",
                                    userSetGpuIds ? "user-" : "auto-",
                                    bPluralGpus ? "s" : "",
                                    numPpRanks,
                                    (numPpRanks > 1) ? "s" : "",
                                    gpuIdsString.c_str());
    }

    /* NOTE: this print is only for and on one physical node */
    GMX_LOG(mdlog.warning).appendText(output);
}

}   // namespace

GpuTaskAssignments::value_type
runTaskAssignment(const MDLogger            &mdlog,
                  const t_commrec           *cr,
                  INodeTaskAssignerPtr       taskAssigner,
                  ConstArrayRef<ComputeTask> computeTasksOnThisRank)
{
    /* Communicate among ranks on this node to find each task that can
     * be executed on a GPU, on each rank. */
    auto tasksOnRanksOfThisNodeEligibleToUseGpus =
        findAllGpuTasksOnThisNode(computeTasksOnThisRank, taskAssigner->areNonbondedOnGpu(),
                                  cr->nrank_pp_intranode, cr->mpi_comm_physicalnode);

    /* Now we have on each rank all the tasks on all ranks of the
     * node, indexed by rank and containing the task type and the ID
     * of assigned GPU. Next we produce the task assignment, if
     * possible. */
    GpuTaskAssignments taskAssignmentOnRanksOfThisNode;
    try
    {
        taskAssignmentOnRanksOfThisNode =
            taskAssigner->assignGpuTasksToDeviceIds(tasksOnRanksOfThisNodeEligibleToUseGpus);
    }
    catch (const std::exception &ex)
    {
        // TODO This implementation is quite similar to that of
        // processExceptionAsFatalError, but it is unclear how we
        // should involve MPI in the implementation of error handling.
        if (cr->rank_intranode == 0)
        {
            printFatalErrorMessage(stderr, ex);
        }
        MPI_Barrier(cr->mpi_comm_physicalnode);
        gmx_exit_on_fatal_error(ExitType_Abort, 1);
    }

    reportGpuUsage(mdlog, taskAssigner->didUserAssignGpus(), taskAssignmentOnRanksOfThisNode,
                   cr->nrank_pp_intranode, cr->nnodes > 1);
    taskAssigner->logPerformanceHints(mdlog, taskAssignmentOnRanksOfThisNode);

    // The work of the taskAssigner is now complete, so release it.
    taskAssigner.release();

    return taskAssignmentOnRanksOfThisNode[cr->rank_intranode];

    // TODO There is no check that mdrun -nb gpu or mdrun -gpu_id is
    // actually being implemented such that nonbonded tasks are being
    // run on compatible GPUs, on all applicable ranks.
}

} // namespace
