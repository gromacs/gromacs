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
 * Defines helper functionality for task assignment on GPUs
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "nonbondedongpu-impl.h"

#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

template <typename T>
size_t countGpuTasksOnThisNode(const std::vector<std::vector<T> > &gpuTasksOnRanksOfThisNode)
{
    size_t numGpuTasksOnThisNode = 0;
    for (const auto &gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        numGpuTasksOnThisNode += gpuTasksOnRank.size();
    }
    return numGpuTasksOnThisNode;
}

//! Instantiate countGpuTasksOnThisNode for GpuTask.
template size_t countGpuTasksOnThisNode<GpuTask>(const std::vector<std::vector<GpuTask> > &);
//! Instantiate countGpuTasksOnThisNode for GpuTaskMapping.
template size_t countGpuTasksOnThisNode<GpuTaskMapping>(const std::vector<std::vector<GpuTaskMapping> > &);

GpuTaskAssignments
buildTaskAssignment(const GpuTasksOnRanks  &gpuTasksOnRanksOfThisNode,
                    const std::vector<int> &flatGpuDeviceIds)
{
    GMX_RELEASE_ASSERT(!flatGpuDeviceIds.empty(), "Must have filled the flat array of device IDs here");
    GpuTaskAssignments gpuTaskAssignmentOnRanksOfThisNode(gpuTasksOnRanksOfThisNode.size());

    auto               currentDeviceId         = flatGpuDeviceIds.begin();
    auto               gpuTaskAssignmentOnRank = gpuTaskAssignmentOnRanksOfThisNode.begin();
    for (const auto &gpuTasksOnRank : gpuTasksOnRanksOfThisNode)
    {
        for (const auto &gpuTaskType : gpuTasksOnRank)
        {
            GMX_RELEASE_ASSERT(currentDeviceId != flatGpuDeviceIds.end(), "Indexing out of range for device IDs");
            (*gpuTaskAssignmentOnRank).push_back({gpuTaskType, *currentDeviceId});
            currentDeviceId++;
        }
        gpuTaskAssignmentOnRank++;
    }
    return gpuTaskAssignmentOnRanksOfThisNode;
}

std::vector<int>
makeFlatGpuTaskAssignment(size_t                  numGpuTasksOnThisNode,
                          const std::vector<int> &compatibleGpus)
{
    std::vector<int> flatGpuDeviceIds;

    // If there are more GPUs detected than tasks, leave the higher-ranked GPUs unused.
    size_t numRanksToShareEachGpu = 1;
    if (numGpuTasksOnThisNode > compatibleGpus.size())
    {
        if (numGpuTasksOnThisNode % compatibleGpus.size() == 0)
        {
            numRanksToShareEachGpu = numGpuTasksOnThisNode / compatibleGpus.size();
        }
        else
        {
            auto message = formatString("The number of GPU tasks (%lu) on this physical node is not a "
                                        "multiple of the number of compatible GPUs detected (%lu). Select "
                                        "a different number of ranks (perhaps specifying PP or PME type), "
                                        "or use the mdrun -gpu_id option to specify manually how the GPUs"
                                        "should be used.", numGpuTasksOnThisNode, compatibleGpus.size());
            GMX_THROW(InconsistentInputError(message));
        }
    }

    // If there are more GPUs than tasks, we leave GPUs empty, and
    // "share" the used ones with 1 rank.

    for (size_t index = 0; index < numGpuTasksOnThisNode; ++index)
    {
        flatGpuDeviceIds.push_back(compatibleGpus[index / numRanksToShareEachGpu]);
    }
    return flatGpuDeviceIds;
}

} // namespace
