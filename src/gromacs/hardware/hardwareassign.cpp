/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "hardwareassign.h"

#include "config.h"

#include <cstring>

#include <algorithm>
#include <numeric>
#include <string>

#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

namespace gmx
{

/*! \brief Check that all user-selected GPUs are compatible.
 *
 * Given the list of selected GPU device IDs in \c gpu_opt and
 * detected GPUs in \c gpu_info, gives a fatal error unless all
 * selected GPUs are compatible
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
 * \param[in]   gpu_info       GPU information including result of compatibility check.
 * \param[in]   gpu_opt        Mapping of GPU IDs derived from the user, e.g. via mdrun -gpu_id.
 */
static void exitUnlessGpuSelectionIsValid(const gmx_gpu_info_t &gpu_info,
                                          const gmx_gpu_opt_t  &gpu_opt)
{
    int         numIncompatibleGpuIds = 0;
    std::string message
        = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    for (int index = 0; index < gpu_opt.n_dev_use; ++index)
    {
        int gpuId = gpu_opt.dev_use[index];
        if (!isGpuCompatible(gpu_info, gpuId))
        {
            numIncompatibleGpuIds++;
            message += gmx::formatString("    GPU #%d: %s\n",
                                         gpuId,
                                         getGpuCompatibilityDescription(gpu_info, gpuId));
        }
    }

    if (numIncompatibleGpuIds > 0)
    {
        gmx_fatal(FARGS, message.c_str());
    }
}

std::vector<int> getCompatibleGpus(const gmx_gpu_info_t &gpu_info)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<int> compatibleGpus;
    compatibleGpus.reserve(gpu_info.n_dev);
    for (int i = 0; i < gpu_info.n_dev; i++)
    {
        GMX_ASSERT(gpu_info.gpu_dev, "Invalid gpu_info.gpu_dev");
        if (isGpuCompatible(gpu_info, i))
        {
            compatibleGpus.push_back(i);
        }
    }
    return compatibleGpus;
}

static std::vector<int> allgather(const int &input, int numRanks, MPI_Comm communicator)
{
    std::vector<int> result(numRanks);
#if GMX_MPI
    MPI_Allgather(&input, 1, MPI_INT,
                  result.data(), 1, MPI_INT,
                  communicator);
#else
    result[0] = data;
#endif
    return result;
}

static std::vector<NonbondedTask> allgatherv(ConstArrayRef<NonbondedTask> input, ConstArrayRef<int> extentOnEachRank, int numRanks, MPI_Comm communicator)
{
    std::vector<int> displacements(numRanks + 1);
    displacements[0] = 0;
    std::partial_sum(std::begin(extentOnEachRank), std::end(extentOnEachRank), std::begin(displacements) + 1);

    // Now allocate the vector and do the allgatherv
    int totalExtent = displacements.back();

    std::vector<NonbondedTask> result(totalExtent);
#if GMX_MPI
    MPI_Allgatherv(input.data(), input.size(), MPI_INT,
                   result.data(), extentOnEachRank.data(), displacements.data(), MPI_INT,
                   communicator);
#else
    result[0] = input;
#endif
    return result;
}

std::vector< gmx::ConstArrayRef<NonbondedTask> >
findAllNonbondedTasksOnThisNode(ConstArrayRef<NonbondedTask> nonbondedTasksOnThisRank,
                                int numRanksOnThisNode,
                                MPI_Comm communicator)
{
    // Find out how many non-bonded tasks are on each rank on this node.
    int numNonbondedTasksOnThisRank = nonbondedTasksOnThisRank.size();
    auto numNonbondedTasksOnEachRankOfThisNode =
        allgather(numNonbondedTasksOnThisRank, numRanksOnThisNode, communicator);

    /* Collect on each rank of this node a vector describing all
     * non-bonded tasks on this node, in ascending order of rank. This
     * requires a vector allgather. */
    auto nonbondedTasksOnThisNode =
        allgatherv(nonbondedTasksOnThisRank, numNonbondedTasksOnEachRankOfThisNode,
                   numRanksOnThisNode, communicator);

    std::vector< gmx::ConstArrayRef<NonbondedTask> > nonbondedTasksOnRanksOfThisNode;
    auto start = std::begin(nonbondedTasksOnThisNode);
    for (int rankId = 0; rankId < numRanksOnThisNode; ++rankId)
    {
        auto end = start + numNonbondedTasksOnEachRankOfThisNode[rankId];
        nonbondedTasksOnRanksOfThisNode.push_back(gmx::ConstArrayRef<NonbondedTask>::fromVector(start, end));
        start = end;
    }
    return nonbondedTasksOnRanksOfThisNode;
}

std::vector< std::vector<NonbondedTask> >
findAllGpuTasksOnThisNode(const std::vector< gmx::ConstArrayRef<NonbondedTask> > &nonbondedTasksOnRanksOfThisNode,
                          int numRanksOnThisNode)
{
    /* Filter the vector of non-bonded tasks on this node to become
     * the vector of non-bonded tasks on ranks of this node that can
     * run on GPUs. */
    std::vector< std::vector<NonbondedTask> > gpuTasksOnRanksOfThisNode(numRanksOnThisNode);
    for (int rankId = 0; rankId < numRanksOnThisNode; ++rankId)
    {
        for (auto taskOnRank : nonbondedTasksOnRanksOfThisNode[rankId])
        {
            // Can this task run on a GPU?
            if (taskOnRank == NonbondedTask::ShortRanged)
            {
                gpuTasksOnRanksOfThisNode[rankId].push_back(taskOnRank);
            }
        }
    }
    return gpuTasksOnRanksOfThisNode;
}

std::vector< std::vector<int> >
mapGpuTasksToDeviceIds(const std::vector< std::vector<NonbondedTask> > &gpuTasksOnRanksOfThisNode,
                       int numRanksOnThisNode,
                       bool forceUsePhysicalGpu,
                       bool tryUsePhysicalGpu,
                       bool rankHasShortRangedTask,
                       bool userSetGpuIds,
                       const gmx_gpu_info_t &gpu_info,
                       const gmx_gpu_opt_t &gpu_opt)
{
    int numGpuTasksOnThisNode = 0;
    for (int rankId = 0; rankId < numRanksOnThisNode; ++rankId)
    {
        numGpuTasksOnThisNode += gpuTasksOnRanksOfThisNode[rankId].size();
    }

    /* Assign IDs of GPUs on this node to the GPU tasks on this
     * node. */
    std::vector< std::vector<int> > gpuDeviceIdsOnRanksOfThisNode(numRanksOnThisNode);
    if (userSetGpuIds)
    {
        if (numGpuTasksOnThisNode != gpu_opt.n_dev_use)
        {
            gmx_fatal(FARGS, gmx::formatString
                      ("You selected %d GPU IDs for the GPU tasks on node %s "
                       "but %lu GPU tasks were found on the ranks of this node",
                       gpu_opt.n_dev_use, getMpiHostname().c_str(), numGpuTasksOnThisNode).c_str());
        }
        // TODO
        /* Validate mdrun -gpu_id, now that it is known from the task
         * assignment whether this rank needs a GPU assigned, which is the
         * user's job if they set -gpu_id. */
        exitUnlessGpuSelectionIsValid(gpu_info, gpu_opt);

        int index = 0;
        for (int rankId = 0; rankId < numRanksOnThisNode; ++rankId)
        {
            for (auto gpuTaskOnRank : gpuTasksOnRanksOfThisNode[rankId])
            {
                GMX_ASSERT(gpuTaskOnRank == NonbondedTask::ShortRanged, "Long-ranged tasks on GPUs are not implemented yet");
                gpuDeviceIdsOnRanksOfThisNode[rankId].push_back(gpu_opt.dev_use[index]);
                index++;
            }
        }
    }
    else
    {
        auto compatibleGpus = getCompatibleGpus(gpu_info);
        if (forceUsePhysicalGpu && compatibleGpus.empty())
        {
            // This logic needs updating when long-ranged tasks can run
            // on GPUs.
            if (rankHasShortRangedTask)
            {
                gmx_fatal(FARGS, "Short-ranged work on GPUs was requested on a rank on node %s, but no compatible "
                          "GPUs were detected. All nodes with such ranks need to have GPUs. If you intended to use "
                          "GPU acceleration in a parallel run, you can either avoid using the nodes that don't have "
                          "GPUs or place PME-only ranks on these nodes.", getMpiHostname().c_str());
            }
        }
        if (tryUsePhysicalGpu && !compatibleGpus.empty())
        {
            if (gpuTasksOnRanksOfThisNode.size() % compatibleGpus.size() != 0)
            {
                gmx_fatal(FARGS, "The number of PP ranks with short-ranged GPU tasks (%lu) on node %s is not "
                          "a multiple of the number of compatible GPUs (%lu). Select a different number of "
                          "ranks, or use the -gpu_id option to manually specify the GPUs to be used.",
                          gpuTasksOnRanksOfThisNode.size(), compatibleGpus.size());
            }
        }

        /* If there are limits on whether ranks can share GPUs,
         * they should be implemented here. */
        auto numRanksToShareEachGpu = gpuTasksOnRanksOfThisNode.size() / compatibleGpus.size();

        int index = 0;
        for (int rankId = 0; rankId < numRanksOnThisNode; ++rankId)
        {
            for (auto gpuTaskOnRank : gpuTasksOnRanksOfThisNode[rankId])
            {
                GMX_ASSERT(gpuTaskOnRank == NonbondedTask::ShortRanged, "Long-ranged tasks on GPUs are not implemented yet");
                gpuDeviceIdsOnRanksOfThisNode[rankId].push_back(index / numRanksToShareEachGpu);
                index++;
            }
        }
    }
    return gpuDeviceIdsOnRanksOfThisNode;
}

DeviceIdForTask getDeviceIdsForThisRank(bool                         rankCanUseGpu,
                                        int numRanksOnThisNode,
                                        bool rankHasShortRangedTask,
                                        const std::vector< std::vector<int> > &deviceIdsOnRanksOfThisNode,
                                        const t_commrec      *cr)
{
    DeviceIdForTask shortRangedTaskGpuIds;

    if (!rankCanUseGpu)
    {
        return shortRangedTaskGpuIds;
    }

    /* Assign IDs of GPUs on this node to the GPU tasks on this
     * node. */
    int index = 0;
    for (int rankId = 0; rankId < numRanksOnThisNode; ++rankId)
    {
        for (auto deviceIdsOnRank : deviceIdsOnRanksOfThisNode[rankId])
        {
            /* Later, there will need to be logic to put the
               device handles into an object whose name or type
               matches what task they will perform. */
            if (rankId == cr->rank_intranode)
            {
                shortRangedTaskGpuIds.push_back(deviceIdsOnRank);
            }
            index++;
        }
    }

    /* Check for a sane assignment given the current capabilities of the code */
    if (rankHasShortRangedTask)
    {
        /* Later, a rank might also have a long-ranged task that could
         * have a GPU handle assigned. */
        GMX_RELEASE_ASSERT(shortRangedTaskGpuIds.size() == 1,
                           gmx::formatString("PP rank %d can use only one GPU, but %lu were assigned",
                                             cr->sim_nodeid, shortRangedTaskGpuIds.size()).c_str());
    }

    return shortRangedTaskGpuIds;
}

} // namespace
