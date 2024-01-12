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
 * Defines routine for collecting all GPU tasks found on ranks of a node.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "findallgputasks.h"

#include "config.h"

#include <numeric>
#include <vector>

#include "gromacs/taskassignment/decidegpuusage.h"
#include "gromacs/taskassignment/taskassignment.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/physicalnodecommunicator.h"

namespace gmx
{

std::vector<GpuTask> findGpuTasksOnThisRank(const bool       haveGpusOnThisPhysicalNode,
                                            const TaskTarget nonbondedTarget,
                                            const TaskTarget pmeTarget,
                                            const TaskTarget bondedTarget,
                                            const TaskTarget updateTarget,
                                            const bool       useGpuForNonbonded,
                                            const bool       useGpuForPme,
                                            const bool       rankHasPpTask,
                                            const bool       rankHasPmeTask)
{
    std::vector<GpuTask> gpuTasksOnThisRank;
    if (rankHasPpTask)
    {
        if (useGpuForNonbonded)
        {
            // Note that any bonded tasks on a GPU always accompany a
            // non-bonded task.
            if (haveGpusOnThisPhysicalNode)
            {
                gpuTasksOnThisRank.push_back(GpuTask::Nonbonded);
            }
            else if (nonbondedTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS,
                          "Cannot run short-ranged nonbonded interactions on a GPU because no GPU "
                          "is detected.");
            }
            else if (bondedTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS,
                          "Cannot run bonded interactions on a GPU because no GPU is detected.");
            }
            else if (updateTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS,
                          "Cannot run coordinate update on a GPU because no GPU is detected.");
            }
        }
    }
    if (rankHasPmeTask)
    {
        if (useGpuForPme)
        {
            if (haveGpusOnThisPhysicalNode)
            {
                gpuTasksOnThisRank.push_back(GpuTask::Pme);
            }
            else if (pmeTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS, "Cannot run PME on a GPU because no GPU is detected.");
            }
        }
    }
    return gpuTasksOnThisRank;
}

namespace
{

//! Constant used to help minimize preprocessing of code.
constexpr bool g_usingMpi = GMX_MPI;

//! Helper function to prepare to all-gather the vector of non-bonded tasks on this node.
std::vector<int> allgather(const int& input, int numRanks, MPI_Comm communicator)
{
    std::vector<int> result(numRanks);
    if (g_usingMpi && numRanks > 1)
    {
        // TODO This works as an MPI_Allgather, but thread-MPI does
        // not implement that. It's only intra-node communication, and
        // happens rarely, so not worth optimizing (yet). Also
        // thread-MPI segfaults with 1 rank.
#if GMX_MPI
        int root = 0;
        // Calling a C API with the const T * from data() doesn't seem
        // to compile warning-free with all versions of MPI headers.
        //
        // TODO Make an allgather template to deal with this nonsense.
        MPI_Gather(const_cast<int*>(&input), 1, MPI_INT, const_cast<int*>(result.data()), 1, MPI_INT, root, communicator);
        MPI_Bcast(const_cast<int*>(result.data()), result.size(), MPI_INT, root, communicator);
#else
        GMX_UNUSED_VALUE(communicator);
#endif
    }
    else
    {
        result[0] = input;
    }

    return result;
}

//! Helper function to compute allgatherv displacements.
std::vector<int> computeDisplacements(ArrayRef<const int> extentOnEachRank, int numRanks)
{
    std::vector<int> displacements(numRanks + 1);
    displacements[0] = 0;
    std::partial_sum(
            std::begin(extentOnEachRank), std::end(extentOnEachRank), std::begin(displacements) + 1);
    return displacements;
}

//! Helper function to all-gather the vector of all GPU tasks on ranks of this node.
std::vector<GpuTask> allgatherv(ArrayRef<const GpuTask> input,
                                ArrayRef<const int>     extentOnEachRank,
                                ArrayRef<const int>     displacementForEachRank,
                                MPI_Comm                communicator)
{
    // Now allocate the vector and do the allgatherv
    int totalExtent = displacementForEachRank.back();

    std::vector<GpuTask> result;
    result.reserve(totalExtent);
    if (g_usingMpi && extentOnEachRank.size() > 1 && totalExtent > 0)
    {
        result.resize(totalExtent);
        // TODO This works as an MPI_Allgatherv, but thread-MPI does
        // not implement that. It's only intra-node communication, and
        // happens rarely, so not worth optimizing (yet). Also
        // thread-MPI segfaults with 1 rank and with zero totalExtent.
#if GMX_MPI
        int root = 0;
        MPI_Gatherv(reinterpret_cast<std::underlying_type_t<GpuTask>*>(const_cast<GpuTask*>(input.data())),
                    input.size(),
                    MPI_INT,
                    reinterpret_cast<std::underlying_type_t<GpuTask>*>(result.data()),
                    const_cast<int*>(extentOnEachRank.data()),
                    const_cast<int*>(displacementForEachRank.data()),
                    MPI_INT,
                    root,
                    communicator);
        MPI_Bcast(reinterpret_cast<std::underlying_type_t<GpuTask>*>(result.data()),
                  result.size(),
                  MPI_INT,
                  root,
                  communicator);
#else
        GMX_UNUSED_VALUE(communicator);
#endif
    }
    else
    {
        for (const auto& gpuTask : input)
        {
            result.push_back(gpuTask);
        }
    }
    return result;
}

} // namespace

/*! \brief Returns container of all tasks on all ranks of this node
 * that are eligible for GPU execution.
 *
 * Perform all necessary communication for preparing for task
 * assignment. Separating this aspect makes it possible to unit test
 * the logic of task assignment. */
GpuTasksOnRanks findAllGpuTasksOnThisNode(ArrayRef<const GpuTask>         gpuTasksOnThisRank,
                                          const PhysicalNodeCommunicator& physicalNodeComm)
{
    int      numRanksOnThisNode = physicalNodeComm.size_;
    MPI_Comm communicator       = physicalNodeComm.comm_;
    // Find out how many GPU tasks are on each rank on this node.
    auto numGpuTasksOnEachRankOfThisNode =
            allgather(gpuTasksOnThisRank.size(), numRanksOnThisNode, communicator);

    /* Collect on each rank of this node a vector describing all
     * GPU tasks on this node, in ascending order of rank. This
     * requires a vector allgather. The displacements indicate where
     * the GPU tasks on each rank of this node start and end within
     * the vector. */
    auto displacementsForEachRank =
            computeDisplacements(numGpuTasksOnEachRankOfThisNode, numRanksOnThisNode);
    auto gpuTasksOnThisNode = allgatherv(
            gpuTasksOnThisRank, numGpuTasksOnEachRankOfThisNode, displacementsForEachRank, communicator);

    /* Next, we re-use the displacements to break up the vector
     * of GPU tasks into something that can be indexed like
     * gpuTasks[rankIndex][taskIndex]. */
    GpuTasksOnRanks gpuTasksOnRanksOfThisNode;
    // TODO This would be nicer if we had a good abstraction for "pair
    // of iterators that point to adjacent container elements" or
    // "iterator that points to the first of a pair of valid adjacent
    // container elements, or end".
    GMX_ASSERT(displacementsForEachRank.size() > 1,
               "Even with one rank, there's always both a start and end displacement");
    auto currentDisplacementIt = displacementsForEachRank.begin();
    auto nextDisplacementIt    = currentDisplacementIt + 1;
    do
    {
        gpuTasksOnRanksOfThisNode.emplace_back(std::vector<GpuTask>());
        for (auto taskOnThisRankIndex = *currentDisplacementIt; taskOnThisRankIndex != *nextDisplacementIt;
             ++taskOnThisRankIndex)
        {
            gpuTasksOnRanksOfThisNode.back().push_back(gpuTasksOnThisNode[taskOnThisRankIndex]);
        }

        currentDisplacementIt = nextDisplacementIt;
        ++nextDisplacementIt;
    } while (nextDisplacementIt != displacementsForEachRank.end());

    return gpuTasksOnRanksOfThisNode;
}

} // namespace gmx
